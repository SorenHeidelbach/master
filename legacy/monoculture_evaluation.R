#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
source("./code/functions.R")
pacman::p_load(
  "gt"
)
checkm <- read_checkm(path = "./mock_ND/monoculture/checkm/storage/bin_stats_ext.tsv") 
checkm2 <- read_checkm(path = "./mock_ND/monoculture/checkm/storage/bin_stats.analyze.tsv") 

checkm  %>% 
  select(1:29) %>% 
  select(bin, marker_lineage, Genome_size, Completeness, Contamination, `0`, `1`, `2`, `3`, `4`, `5+`, GC) %>% 
  mutate(
    bin = str_remove(bin, "MinION_"),
    marker_lineage = marker_lineage,
    GC = signif(as.numeric(GC), 2)*100,
    Contamination = signif(as.numeric(Contamination), 2),
    Completeness = signif(as.numeric(Completeness), 2),
    Genome_size = scales::scientific(as.numeric(Genome_size))
  ) %>% 
  rename(c("Sample"="bin")) %>% 
  arrange(Sample) %>% 
  gt()  %>%
  data_color(
    columns = Completeness,
    colors = scales::col_numeric(
      palette = c(
        "red1", "orange",  "yellow1", "green1"),
      domain = c(20, 100))
  ) %>%
  data_color(
    columns = Contamination,
    colors = scales::col_numeric(
      palette = c(
        "green1", "yellow1", "red1"),
      domain = c(4, 0))
  )

read_gfa <- function(gfa_path){
  system(paste0("awk '/^P/ {print $2,$3}' ", gfa_path), intern = TRUE ) %>% 
    data.frame() %>% 
    separate(col = ".", into = c("contig", "edge"), sep = " ") %>% 
    mutate(
      edge = str_remove(edge, "\\+")
    ) %>% 
    left_join(
      system(paste0("awk '/^S/ {print $2,$4}' ", gfa_path), intern = TRUE ) %>% 
        data.frame() %>% 
        separate(col = ".", into = c("edge", "cov"), sep = " ") %>% 
        mutate(cov = str_remove(cov, "dp\\:i\\:"))
    )
}

read_gfa("/shared-nfs/SH/mock_ND/monoculture/flye/MinION_BA_NAT/assembly_graph.gfa")
system("echo /shared-nfs/SH/mock_ND/monoculture/flye/*", intern = TRUE) %>% 
  str_split(" ") %>% 
  unlist() %>% 
  paste0("/assembly_graph.gfa") %>% 
  map(read_gfa) %>% view

