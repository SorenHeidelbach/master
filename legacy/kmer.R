#!/usr/bin/Rscript

# Setting up library path and laoding functions
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
source("./code/functions.R")
pacman::p_load(
  "tidyverse",
  "ape",
  "keras"
)
# Creating directionary for figures
dir.create("./figures/kmer_clustering")

# Loading nanodisco motifs
MGM1_motif <- load_ND_profile("/shared-nfs/SH/nanodisco/MGM1/dataset/methylation_profile_MGM1_motif.RDS")
MGM1_motif_all <- readRDS("/shared-nfs/SH/nanodisco/MGM1/dataset/")
# Loading coverage dataframe and calculating GC%
coverage <- attributes(readRDS("./nanodisco/MGM1/dataset/methylation_profile_MGM1_motif.RDS"))$contig_coverage %>% 
  left_join(get_GC_percent("./nanodisco/MGM1/reference/metagenome.fasta"), by = c("chr" = "contig")) %>% 
  filter(chr %in% (MGM1_motif$contig %>% unique))

# Converting motifs dataframe into kmer features
MGM1_feature <- MGM1_motif %>% 
  unnest(cols = c(data)) %>%
  group_by(contig, motif) %>% 
  dplyr::summarise(
    nb_occurrence = mean(nb_occurrence)
  ) %>%
  select(contig, motif, nb_occurrence) %>%
  ungroup() %>% 
  pivot_wider(values_from = nb_occurrence, names_from = motif) %>% 
  mutate_if(
    is.numeric,
    function(x) x = ifelse(is.na(x), 0, as.numeric(x))
    ) 

##-----------------------------------------------------------------------------
# kmer calculation
metagenome_fasta <- seqinr::read.fasta("./nanodisco/MGM1/reference/metagenome.fasta" )
data.frame(
  contig_names = (metagenome_fasta %>% 
    `[`(metagenome_fasta %>% map_dbl(length) > 25000) %>% 
    names())
) %>% 
  mutate(
    seq = metagenome_fasta[contig_names]
  )
kmer_length = 8
kmer_freq <- metagenome_fasta %>% 
  `[`(metagenome_fasta %>% map_dbl(length) > 25000) %>% 
  names() %>% 
  map(function(contig) {
    print(contig)
    seq <- metagenome_fasta[contig] %>% 
      unlist() %>% 
      paste0(collapse = "")

    1:(nchar(seq)-kmer_length) %>% 
      map(function(x) {
        substr(seq, start = x, stop = (x+kmer_length))
      }) %>% 
      unlist %>% 
      table %>% 
      data.frame %>% 
      set_names(c("kmer", contig))
    }) %>% 
  reduce(left_join, by = "kmer") %>% 
  mutate_if(is.numeric,
    function(x) ifelse(is.na(x), 0, x)
  ) %>% 
  mutate_if(is.numeric, 
    function(x) (x - min(x)) / (max(x) - min(x))          
  ) %>% 
  `rownames<-`(.$kmer)
test <- kmer_freq %>% 
  select(-kmer) %>% 
  data.table::as.data.table() %>% 
  data.table::transpose() %>% 
  set_names(kmer_freq$kmer) %>% 
  as.data.frame() %>% 
  tsne() 
ggplot(test, aes(x = V1, y = V2)) +
  geom_point()

##-----------------------------------------------------------------------------
# PCA

# PC transformation
MGM1_feature_pca <- MGM1_feature %>%  
  select(-contig) %>% 
  pca() %>% 
  mutate(
    cluster = hdbscan(cbind(PC1, PC2), minPts = 13)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = MGM1_feature$contig) %>% 
  left_join(coverage)
  
# Scatter plot with contigs
ggplot(MGM1_feature_pca, aes(x = PC1, y = PC2, size = contig_length.x)) +
  theme_bw() +
  stat_density_2d(col = "gray60") +
  geom_point(aes(fill = cluster), shape = 21, alpha = 0.6) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA"
  )

##-----------------------------------------------------------------------------
# t-SNE

# t-SNE transformation 
MGM1_feature_tsne <- MGM1_feature %>%  
  select(-contig) %>% 
  standardise_cols() %>% 
  tsne() %>% 
  mutate(
    cluster = hdbscan(cbind(V1, V2), minPts = 13)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = MGM1_feature$contig) %>% 
  left_join(coverage)

# Scatter plot with contigs
ggplot(MGM1_feature_tsne, aes(x = V1, y = V2, size = contig_length.x)) +
  theme_bw() +
  stat_density_2d(col = "gray60") +
  geom_point(aes(fill = cluster), shape = 21, alpha = 0.6) +
  labs(
    x = "dim 1",
    y = "dim 2",
    title = "t-SNE"
  )

##-----------------------------------------------------------------------------
# Autoencoder

# Auto encoder latent variable extraction
MGM1_feature_ae <- MGM1_feature %>%  
  select(-contig) %>% 
  normalise_cols() %>% 
  autoencoder() %>% 
  mutate(
    cluster = hdbscan(cbind(V1, V2), minPts = 13)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = MGM1_feature$contig) %>% 
  left_join(coverage)

# Scatter plot with contigs
ggplot(MGM1_feature_ae, aes(x = V1, y = V2, size = contig_length.x)) +
  theme_bw() +
  stat_density_2d(col = "gray60") +
  geom_point(aes(fill = cluster), shape = 21, alpha = 0.6) +
  labs(
    x = "dim 1",
    y = "dim 2",
    title = ""
  )
