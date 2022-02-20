

.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
pacman::p_load(
  "data.table",
  "dplyr",
  "ggplot2",
  "rhdf5",
  "ggforce"
)

read <- H5Fopen("samples/zymo/WGA/fast5_single_filt/0/000aaa32-5a93-41b2-b7b8-4e6f6cff3d0c.fast5")

h5ls(read)

read$"/Analyses/RawGenomeCorrected_002/BaseCalled_template/Events" %>% 
  ggplot(aes(x = start, y = norm_mean, col = base)) +
  geom_line(aes(group = 1)) +
  facet_zoom(xlim = c(1e4, 1.05e4))
