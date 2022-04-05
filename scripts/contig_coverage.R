pacman::p_load(
  "ggplot2",
  "dplyr",
  "tidyr",
  "data.table",
  "stringr",
  "argparse",
  "seqinr"
)
###############################################################################
# Arguments
###############################################################################


arg <- list()
# File generated from 'SAMtools depth -a' 
arg$cov_file <- "/scratch/users/sheide17/samples/zymoHMW/PCR/trim_filt_ref.depth"
arg$reference <- "/shared-nfs/SH/samples/zymoHMW/zymoHMW_references.fasta"
arg$out <- "/shared-nfs/SH/samples/zymoHMW/pcr"

dir.create(arg$out)

reference <- read.fasta(arg$reference) %>% 
  lapply(function(x){
    list(
      attributes(x)$name,
      length(x),
      sum(grepl("[gGcC]", x))/length(x)
      )
  }) %>% 
  rbindlist() %>% 
  setnames(c("ref", "size", "GC"))


n = 50

coverage <- fread(arg$cov_file) %>% 
  setnames(c("ref", "pos", "cov"))
coverage_binned <- coverage[
      , bin := c(rep(1:(floor(.N/n)), each = n), rep(floor(.N/n)+1, .N %% n)), by = ref
    ][
      , .(cov_bin = as.integer(mean(cov)), start_pos = min(pos)), by = .(bin, ref)
    ][
      reference, on = .(ref), `:=`(size = i.size, GC = i.GC)
    ][
      , ref := str_replace(ref, "E[a-z]*_", "E\\."),
    ][
      , ref := str_remove(ref, "_genome"),
    ][
      , ref := str_remove(ref, "_chromosome"),
    ][
      , ref := str_remove(ref, "_complete"),
    ][
      , ref := str_replace(ref, "_plasmid", " p"),
    ][
      , ref := str_remove(ref, "\\.pilon.*"),
    ][
      , ref := str_replace(ref, "S[a-z]*_", "S\\."),
    ][
      , ref := str_replace(ref, "L[a-z]*_", "L\\."),
    ][
      , ref := str_replace(ref, "P[a-z]*_", "P\\.")
    ][
      , ref := str_replace(ref, "tig000.*", "S.cerevisiae")
    ][
      , ref := str_replace(ref, "BS", "B.subtilis")
    ][
      , GC := mean(unique(GC)), by = ref
    ]


coverage_binned %>% 
  mutate(cov_bin = ifelse(cov_bin == 0, cov_bin + 0.5, cov_bin)) %>% 
  ggplot() +
  aes(y = cov_bin, x = ref) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_text(
    data = coverage_binned[, .(median = median(cov_bin)), by = ref], 
    aes(x = ref, y = median, label = median), hjust = 1) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 0)) +
  scale_y_log10() +
  labs(
    x = "Reference",
    y = "Coverage"
  )

ggsave(
  filename = paste0(arg$out, "/coverage.png"),
  device = "png",
  width = unit(9, "cm"),
  height = unit(6, "cm"),
  coverage_binned %>% 
    mutate(cov_bin = ifelse(cov_bin == 0, cov_bin + 0.5, cov_bin)) %>% 
    ggplot() +
    aes(x = cov_bin, y = ref) +
    geom_violin(aes(fill = GC), scale = "width", orientation = "y", adjust = 1.5) +
    geom_text(
      data = coverage_binned[, .(median = median(cov_bin)), by = ref], 
      aes(y = ref, x = median, label = median), hjust = 0.5) +
    theme_bw() +
    theme(axis.text.y = element_text(hjust = 1)) +
    scale_x_log10() +
    labs(
      y = "Reference",
      x = "Coverage"
  ) +
    scale_fill_gradient2(low = "#F7977A", mid = "#FFF79A", high = "#82CA9D", midpoint = 0.5)
)

ggsave(
  filename = paste0(arg$out, "/proportion_mapped.png"),
  device = "png",
  width = unit(7, "cm"),
  height = unit(6, "cm"),
  coverage[
    , .(bases = sum(cov)), by = ref
  ][
    , proportion := bases / sum(bases)
  ] %>% 
    filter(!grepl("tig000", ref)) %>% 
    ggplot(aes(x = ref, y = proportion)) +
    geom_histogram(stat = "identity") +
    theme_bw() +
    theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1)) +
    geom_text(aes(label = paste0(signif(bases/1e6, 3), "Mb")), nudge_y = 0.005) +
    labs(
      x = "Reference",
      y = "Proportion of bases mapped to reference"
    )
)

coverage_binned %>% 
  ggplot() +
  aes(x = bin*n, y = cov_bin) +
  geom_line() +
  facet_wrap(~ref, scales = "free") +
  theme_bw() +
  labs(x = "Contig position", y = "Bin coverage", title = paste0("Reference coverage, bin size ", n)) +
  geom_smooth(col = "red3")



