
# Setting costum library path
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
if (!require("pacman")) install.packages("pacman", lib = "/shared-nfs/SH/code/Rlib")
if (!require("BiocManager")) install.packages("BiocManager", lib = "/shared-nfs/SH/code/Rlib")

# Loading libraries with pacman
pacman::p_load( 
  "dplyr",
  "ggplot2",
  "ggforce",
  "data.table"
)

###############################################################################
## Preprocessing of Nanopolish Event Align
###############################################################################

event <- fread("/shared-nfs/SH/samples/zymo/WGA_event_align.txt")

event_mean <- event[
    , .(event_level_mean_pos = mean(model_mean - event_level_mean)), by = .(contig, position, reference_kmer)
  ][
    , NT := strtrim(reference_kmer, width = 1)
  ]
ggplot(pos_mean, aes( x = position, y = event_level_mean_pos, col = NT)) +
  geom_line(size = 1, alpha = 0.5) +
  facet_zoom(x = c(50000, 50100)) +
  theme_bw()

values_relative_position <- function(DT, query_motif, frame = -10:10) {
  DT[
    , .(event_level_mean_pos = mean(event_level_mean)), by = .(contig, position, reference_kmer)
    ][
    , motif := grepl(paste0("^", query_motif), reference_kmer)
    ][
    , paste0("p_", frame) := lapply(frame, shift, x = event_level_mean_pos), by = contig
    ][
      (motif)
    ][
      , motif_seq := query_motif
    ] %>% 
      melt.data.table(measure.vars = patterns("^p_"), 
                      variable.name = "rel_pos")
}
event_motif <- lapply(c("AAAAAA", "TTTTTT"), values_relative_position, DT = event) %>% 
  rbindlist() 

ggplot(event_motif, aes(x=rel_pos, y = value)) +
  geom_jitter(alpha = 0.05, size = 0.5) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  facet_wrap(~motif_seq)

