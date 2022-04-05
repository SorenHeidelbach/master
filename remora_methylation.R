.libPaths(c("/shared-nfs/SH/scripts/Rlib", .libPaths()))
Sys.setenv("R_LIBS_USER" = "Rlib")
setwd(here::here())

pacman::p_load(
  "data.table",
  "ggplot2",
  "dplyr",
  "stringr",
  "seqinr"
  )

# Reference GC %
arg <-  list()
arg$reference <- "/shared-nfs/SH/data/zymoHMW/zymoHMW_references.fasta"
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


# Load modbed file from remora
remora <- fread("/shared-nfs/SH/data/zymoHMW/NAT/bonito/basecalls_v2_sort2.bed") %>% 
  setnames(c("ref", "chr_start", "chr_end", "mod_name", "score", "strand", "drop1", "drop2", "drop3", "coverage", "percent_read_mod"))

remora <- remora[
    reference, on = .(ref), `:=`(size = i.size, GC = i.GC)
  ][
  , .SD, .SDcols = ! names(remora) %like% "drop"
  ][
  , ref := str_replace(ref, "E[a-z]*_", "E\\. "),
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
  , ref := str_replace(ref, "S[a-z]*_", "S\\. "),
  ][
  , ref := str_replace(ref, "L[a-z]*_", "L\\. "),
  ][
  , ref := str_replace(ref, "P[a-z]*_", "P\\. ")
  ][
  , ref := str_replace(ref, "tig000.*", "S. cerevisiae")
  ][
  , ref := str_replace(ref, "BS", "B. subtilis")
  ][
  , ref := paste0("*", ref, "*")
  ]

lapply(
  seq(25, 100, by = 5),
  function(x){
    remora[
        percent_read_mod > x & ! grepl(pattern = "cerevisiae", x = ref) & ! grepl(pattern = " p", x = ref)  , .(count = .N, rel_count = .N/(max(GC)*max(size)*2)), by = .(ref, strand)
      ][
        , threshold := x
      ]
  }
) %>%  rbindlist() %>% 
  ggplot() +
  aes(x = threshold, y = rel_count, col = ref) +
  geom_line() +
  facet_wrap(~strand) +
  theme_bw() +
  labs(
    x = "Percent Mapped Reads with Mod",
    y = "Percent of C's with mod"
  )
remora[
    percent_read_mod > 80, .(count = .N), by = ref
  ] 

remora[
  ref == "*E. coli*" & coverage > 10 & chr_start > 10 & chr_start < 50
] %>% 
  ggplot() +
  aes(x = chr_start, percent_read_mod) +
  geom_point(size = 0.5) +
  theme_bw()










