#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
pacman::p_load(
  "dplyr",
  "stringr",
  "rhdf5",
  "purrr"
)
arg <- commandArgs()
path <- arg[6]
out <- arg[7] %>% 
  paste0(
    "/",
    basename(path) %>% 
      str_remove(".fast5") %>% 
      paste0(".fastq")
  )

h5ls(path) %>%
  filter(str_detect(group, "BaseCalled_template$") & name == "Fastq") %>%
  pull(group) %>%
  map_chr(function(x){
    h5read(path, name = x)$Fastq
  }) %>%
  write(out, append = TRUE)


