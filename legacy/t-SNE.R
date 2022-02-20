
# t-SNE feature transformation
tsne_features <- 
  features %>% 
  select(-contig) %>% 
  standardise_cols() %>% 
  tsne() %>% 
  mutate(
    cluster = hdbscan(cbind(V1, V2), minPts = p$hdbscan$minPts)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% 
      factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = features$contig) %>% 
  left_join(coverage) %>% 
  group_by(cluster) %>% 
  mutate(
    bin_size = sum(contig_length.x),
    n_contigs = length(cluster)
  )
write_csv(tsne_features, file = paste0(p$misc$output_folder, "/features_tsne.csv"))

# Plot tSNE with clustering
ggsave(
  filename = paste0(p$misc$output_folder, "/cluster_tsne.pdf"),  
  width = 8, height = 6,
  ggplot(
    tsne_features, 
    aes(x=V1, y=V2, fill=cluster, size = contig_length.x)
  ) +
    geom_point(
      aes(fill = cluster), 
      shape = 21, 
      alpha = 0.8
    ) +
    scale_size(
      range = c(1, 15), 
      limits = c(25000, 900000)
    ) +
    scale_fill_manual(
      values = c(
      '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
      '#fabed4', '#469990', '#dcbeff', '#fffac8', 'gray50', 
      '#aaffc3', '#000075', '#a9a9a9'), na.value = "white"
    ) +
    labs(
      x = "t-SNE dimension 1", 
      y = "t-SNE dimension 2", 
      title = "Contig Clustering"
    ) +
    guides(
      size = guide_legend(title = "Contig Size"),
      fill = guide_legend(title = "Cluster", override.aes = list(size = 5))
    ) +
    theme_bw()
)

##-----------------------------------------------------------------------------
# GC vs. cov but with all cluster added
##-----------------------------------------------------------------------------
ggsave(
  filename = paste0(p$misc$output_folder, "/gc_vs_cov_tsne_clusters.pdf"),  
  width = 8, height = 6,
  ggplot(tsne_features, aes(x=avg_cov.dataset_B, y=gc, fill=cluster)) +
    geom_point(aes(size = contig_length.x), alpha = 0.4, shape = 21) +
    scale_size(range = c(1, 15), limits = c(25000, 900000), 
               guide = guide_legend(title = "Contig Length")) +
    scale_x_log10() +
    labs(y = "GC%", x = "Mean Contig Coverage") +
    theme_bw()
)

# Saving cluster fasta
write_cluster_fasta(
  df = tsne_features, 
  metagenome_fasta = fasta,
  path_out = paste0(p$misc$output_folder, "/bin_fasta_tSNE/")
)
