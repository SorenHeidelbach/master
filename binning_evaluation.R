
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
pacman::p_load(
  "data.table",
  "tidyr"
)


metabat <- fread("samples/zymo/metabat2/test.tsv")
