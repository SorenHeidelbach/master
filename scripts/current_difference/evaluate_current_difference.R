#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/software/Rlib" , .libPaths()))
pacman::p_load(
  "data.table",
  "dplyr",
  "ggplot2", 
  "seqinr",
  "stringr",
  "patchwork",
  "logger",
  "glue",
  "pbapply",
  "scales",
  "ggforce",
  "argparse",
  "tidyr"
)
setwd(here::here())

check_adjacancy_vec <- function(vec1, vec2, min_dist = 10){
  # Check if any values in vec1 are within 'min_dist' of any values in vec2
  dist(c(vec2,vec1), upper = TRUE) %>% 
    as.matrix() %>% 
    `[`((length(vec2)+1):(length(vec2)+length(vec1)), 1:length(vec2)) %>% 
    apply(MARGIN = 1, function(z) any(z <= min_dist))
}

generate_kmers <- function(kmer_size, capital = FALSE){
  # Generate vector of all possible kmers of size 'kmer_size'
  i <- 2
  kmers <- c("a", "g", "c", "t")
  while (i <= kmer_size) {
    kmers <- expand.grid(kmers, c("a", "g", "c", "t")) %>% apply(MARGIN = 1, function(n) paste0(n, collapse = ""))
    i <- i +1
  }
  if ( capital){
    kmers <- toupper(kmers)
  }
  return(kmers)
}

# Arguments
###############################################################################
# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--data_dir", help='Directionary containing output from "calculate_current_difference.R"', required = TRUE)
parser$add_argument("--reference",  nargs=1, help="Reference mapped too", required = TRUE)
parser$add_argument("--out",  nargs=1, help="Output folder [default data_dir]")
arg <- parser$parse_args()

arg$min_p_val <- 1e-12
# arg$data_dir <- "/shared-nfs/SH/results/zymoHMW/current_difference/test1_20220331"
# arg$out <- "/shared-nfs/SH/results/zymoHMW/current_difference/test1_20220331"
if (is.null(arg$out)) {
  arg$out <- arg$data_dir
}

arg <-  list()
arg$reference <- "/shared-nfs/SH/data/zymoHMW/zymoHMW_references.fasta"
arg$data_dir <- glue("/shared-nfs/SH/results/zymoHMW/current_difference/test1_20220404")




current_difference <- fread(glue("{arg$data_dir}/current_difference.tsv.gz"))
current_difference[
  , event := p_val < arg$min_p_val 
  ][
  , distance_previous_event := contig_index - data.table::shift(contig_index, 1), by = .(contig, event, direction)
  ][
  , distance_next_event := data.table::shift(distance_previous_event, -1), by = .(contig, event, direction)
  ][
  , event_group := if_else(distance_previous_event > 2 | is.na(distance_previous_event), seq_len(.N), NULL)
  ][
  , event_group := event_group[1], by = .(cumsum(!is.na(event_group)), event)
  ][
  , event_center := contig_index == contig_index[ceiling(.N/2)], by = .(contig, event_group)
  ][
  , event_center := if_else(event, event_center, event)
  ]
current_difference_pcr <- fread(glue("{arg$data_dir}/current_difference_pcr.tsv.gz"))
current_difference_pcr[
  , event := p_val_weighted_log < arg$min_p_val 
  ][
  , distance_previous_event := contig_index - data.table::shift(contig_index, 1), by = .(contig, event, direction)
  ][
  , distance_next_event := data.table::shift(distance_previous_event, -1), by = .(contig, event, direction)
  ][
  , event_group := if_else(distance_previous_event > 2 | is.na(distance_previous_event), seq_len(.N), NULL)
  ][
  , event_group := event_group[1], by = .(cumsum(!is.na(event_group)), event)
  ][
  , event_center := contig_index == contig_index[ceiling(.N/2)], by = .(contig, event_group)
  ][
  , event_center := if_else(event, event_center, event)
  ]
# Plots 
pboptions(type = "timer")
plot_events <- function(dt){
  guides_format <- guides(
    fill = guide_legend(title = "Event"),
    size = guide_legend(title = "Event"),
    col = guide_legend(title = "Event"),
    alpha = "none"
  )
  plot_col <- c("#99712d", "#2d5599")
  geom_point_default <- geom_jitter(alpha = 0.1, size = 0.1)
  
  
  p1.1 <- dt %>% 
    ggplot(aes(x = contig_index, y = mean_dif, fill = event_center)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = mean_dif), y = 0, size = 0.2) +
    geom_point(shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (Mean difference)"
    )  + 
    scale_size_manual(values = c(0.4, 2)) +
    scale_fill_manual(values = plot_col) +
    guides_format +
    ggtitle("mean_diff")
  
  p1.2 <- dt %>% 
    ggplot(aes(x = contig_index, y = p_val_weighted_log)) +
    geom_hline(yintercept = 1) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = p_val_weighted_log), y = 0, size = 0.2) +
    geom_point(shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (U-value)"
    ) +
    scale_size_manual(values = c(0.4, 2)) +
    scale_y_log10() +
    scale_fill_manual(values = plot_col) +
    guides_format +
    ggtitle("p_val")
  
  
  p1.3 <- dt %>%  
    ggplot(aes(x = n_pcr, y = n_nat)) +
    labs(
      x = "PCR dacs",
      y = "NAT dacs"
    ) + 
    geom_point(size = 0.3, alpha = 0.1, col = "gray30") +
    stat_density_2d_filled(
      contour = TRUE,
      aes(alpha = after_stat(level), fill = after_stat(level)),
      bins = 10, col = "gray10", size = 0.1
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_alpha_discrete(range = c(0,1)) +
    scale_fill_manual( values = c(colorRampPalette(c("red",  "darkred"))(10))) +
    ggtitle("dacs") +
    theme_bw()
  
  p1.4 <- dt %>% 
    ggplot(aes(x = n_pcr_map, y = n_nat_map)) +
    labs(
      x = "PCR coverage",
      y = "NAT coverage"
    ) + 
    geom_jitter(size = 0.3, alpha = 0.1, col = "gray30") +
    stat_density_2d_filled(
      contour = TRUE,
      aes(alpha = after_stat(level), fill = after_stat(level)),
      bins = 10, col = "gray10", size = 0.1
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_alpha_discrete(range = c(0,1)) +
    scale_fill_manual(values = c(colorRampPalette(c("red",  "darkred"))(10))) +
    theme_bw() +
    ggtitle("coverage")
  
  p1.5 <- dt %>% 
    mutate(
      mappings = n_pcr_map / n_nat_map
    ) %>% 
    ggplot(aes(y = mappings, x = p_val_weighted_log)) +
    labs(
      x = "U-value",
      y = "nr. PCR / nr. NAT"
    ) + 
    geom_point_default +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format +
    ggtitle("cov_proportion")
  
  p1.6 <- dt %>% 
    ggplot(aes(x = mean_dif, y = p_val_weighted_log)) +
    labs(
      x = "Mean difference",
      y = "U-value"
    ) + 
    geom_point_default +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format +
    ggtitle("fake_vulcano_plot")
  
  return(list(p1.1, p1.2, p1.3, p1.4, p1.5, p1.6))
}

p_all <- plot_events(current_difference)
p_all %>%
  pblapply(
    function(p){
      ggsave(
        width = unit(8, "cm"),
        height = unit(10, "cm"),
        p,
        filename = glue("{arg$data_dir}/current_difference_{p$labels$title}.png"),
        device = "png"
      )
    }
  )


plot_col <- c("red", "#ffa600")
ggsave(
  width = unit(8, "cm"),
  height = unit(10, "cm"),
  filename = glue("{arg$out}/current_difference_zoom.png"),
  device = "png",
  current_difference %>% 
    mutate(mean_dif2 = abs(dist) * if_else(direction == "fwd", 1, -1)) %>% 
    ggplot() +
    aes(x = contig_index, y = mean_dif2) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = mean_dif2, col = direction, alpha = event_center), y = 0, size = 0.1) +
    #geom_point(aes(size = event_center, fill = event_center), shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (Mean difference)"
    )  + 
    scale_size_manual(values = c(0.2, 2)) +
    scale_fill_manual(values = plot_col) +
    scale_color_manual(values = c("#1034a6", "#850101")) +
    scale_alpha_manual(values = c(0.01, 1)) +
    guides(
      fill = guide_legend(title = "Event"),
      size = guide_legend(title = "Event"),
      col = guide_legend(title = "Direction"),
      alpha = "none"
    ) +
    theme_bw()
)


# plot_events(signal_mappings_pcr_0[
#   n_nat_map > 1 & n_pcr_map > 1
# ])

# ggsave(
#   width = unit(8, "cm"),
#   height = unit(10, "cm"),
#   p_pcr_vs_pcr, 
#   filename = glue("{arg$out}/current_difference_pcr_vs_pcr_{arg$contig_plot}.png"), 
#   device = "png"
# )



###############################################################################
# Evaluation of differenct U-value cut offs
###############################################################################

# pcr vs. pcr and nat vs. pcr plot function
plot_PvP_NvP <- function(thresholds, col){
  thresholds %>% 
    lapply(
      function(x){
        n_NvP <- sum(abs(current_difference[, ..col]) < x, na.rm = TRUE)
        p_NvP <- n_NvP / nrow(na.omit(current_difference[, ..col]))
        n_PvP <- sum(abs(current_difference_pcr[, ..col]) < x, na.rm = TRUE)
        p_PvP <- n_PvP / nrow(na.omit(current_difference_pcr[, ..col]))
        list(
          x,
          p_NvP, 
          p_PvP,
          p_NvP - p_PvP,
          p_PvP/p_NvP
        )
      }) %>% 
    rbindlist() %>% 
    setnames(paste0("V", 1:5), c("P.val", "p_NvP",  "p_PvP", "p_diff", "p_ratio")) %>%  
    tidyr::pivot_longer(starts_with("p_"), values_to = "proportion", names_to = "type") %>% 
    ggplot() +
    aes(x = P.val, y = proportion, col = type) +
    geom_line(size = 1.2) +
    #geom_point(shape = 21, size = 1.5) +
    scale_x_log10() +
    scale_y_continuous(labels = percent) +
    labs(
      x = "P-value threshold",
      y = "Proportion of events"
    ) +
    theme_bw() +
    scale_color_manual(values = c("#127474", "#741212", "#437412", "#431274"))
}
ylimits <- c(1e-5, 1)
ybreaks <- 10 ^ seq(log(ylimits[2], base = 10), log(ylimits[1], base = 10), by = -1)
ylabel <- paste0(signif(ybreaks*100, 1), "%")
# p-value
p2.1 <-  10^-seq(0, 2, by = 0.2) %>% 
  plot_PvP_NvP(col = "p_val") +
  labs(title = "P-value") 
  #scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

# p-value log weighted 
p2.3 <-  10^-seq(0, 2, by = 0.2) %>% 
  plot_PvP_NvP(col = "p_val_weighted_log") +
  labs(title = "P-value log weighted") 
  #scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

# Mean difference cutoff
p2.4 <- 10^-seq(-1, 3, by = 0.1) %>% 
  plot_PvP_NvP("dist") +
  xlab("Mean difference threshold") +
  labs(title = "Mean difference")

p2.5 <- 10^-seq(0, 2, by = 0.1) %>% 
  plot_PvP_NvP("p_val_BH") +
  labs(title = "P-value B&H corrected") 
  #scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

p2.6 <- 10^-seq(0, 2, by = 0.1) %>% 
  plot_PvP_NvP("p_val_BH_weighted_log") +
  labs(title = "P-value B&H corrected and log weighted") 
  #scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

p2.7 <- 10^-seq(0, 2, by = 0.1) %>% 
  plot_PvP_NvP("p_val_bonf") +
  labs(title = "P-value Bonferroni") 
  #scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)


p2_all <- (p2.1 + p2.5 + p2.7) /
  (p2.3 + p2.6 + p2.4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title =  "Estimation of event cut-off"
  ) & 
  theme_bw()
ggsave(
  width = unit(16, "cm"),
  height = unit(7, "cm"),
  filename = glue("{arg$out}/event_cutoff_estimation.png"),
  p2_all,
  device = "png")



###############################################################################
# Motif event enrichment
###############################################################################

reference <- read.fasta(arg$reference)

reference <- reference[["Escherichia_coli_plasmid"]] %>% 
  unlist() %>% 
  paste0(collapse = "") 
event_center <- current_difference[event_center == TRUE, ]$contig_index
event_center_pcr <- current_difference_pcr[event_center == TRUE, ]$contig_index
motifs <- generate_kmers(4)

motif_event_frequency <- lapply(
  motifs,
  function(x){
    str_locate_all(reference, x) %>% 
      as.data.frame %>% 
      mutate(
        motif = x,
        event_adjacent =  check_adjacancy_vec(start, event_center, min_dist = 10),
        event_adjacent_pcr = check_adjacancy_vec(start, event_center_pcr, min_dist = 10)
        )
  }     
  ) %>% 
  rbindlist() %>% 
  group_by(motif) %>% 
  summarise(
    kmer_count = n(),
    percent_event_adjacent = sum(event_adjacent)/kmer_count,
    percent_event_adjacent_pcr = sum(event_adjacent_pcr)/kmer_count
  ) %>% 
  pivot_longer(cols = c(percent_event_adjacent, percent_event_adjacent_pcr), names_to = "type", values_to = "percent") %>% 
  mutate(
    motif = factor(motif, levels = motif[order(percent, decreasing = TRUE)] %>% unique())
  )
y_scale <- max(motif_event_frequency$kmer_count)
ggplot(motif_event_frequency) +
  aes(x = motif, fill = type) +
  geom_line(aes(y = 0.04 * kmer_count/y_scale), group = 1, col = "black") +
  geom_point(aes(y = percent), shape = 21) +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, size = 8)) +
  facet_wrap(~type) +
  scale_y_continuous(
    name = "Percent of kmers adjacent to event",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./0.04*y_scale, name="Number of kmer in reference")
  )


lapply(
  motifs,
  function(x){
    str_locate_all(reference, x) %>% 
      as.data.frame %>% 
      mutate(
        motif = x
      )
  }     
) %>% 
  rbindlist() %>% 
  ggplot() +
  aes(x = motif) +
  geom_histogram(stat = "count")
