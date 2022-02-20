
ae_features <- features %>% 
  select(-contig) %>% 
  autoencoder(
    activator = p$autoencoder$activation_function, 
    intermediate_nodes = p$autoencoder$intermediate_nodes,
    bottleneck_nodes = p$autoencoder$bottleneck_nodes,
    epochs = p$autoencoder$epochs
  )  %>% 
  select(starts_with("V")) %>% 
  mutate(
    cluster = hdbscan(., minPts = p$hdbscan$minPts)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = features$contig) %>% 
  left_join(coverage) %>% 
  group_by(cluster) %>% 
  mutate(
    bin_size = sum(contig_length.x),
    n_contigs = length(cluster)
  )
write_csv(ae_features, file = paste0(p$misc$output_folder, "/features_AE.csv"))

ggsave(
  filename = paste0(p$misc$output_folder, "/cluster_AE.pdf"),  
  width = 8, height = 6,
  ggplot(
    ae_features,
    aes(x = V1, y = V2, size = contig_length.x)) + 
    geom_density_2d(
      contour = TRUE, 
      col = "gray60") +
    geom_point(
      aes(fill = cluster), 
      shape = 21, 
      alpha = 0.6) +
    theme_bw() +
    labs(
      x = "Latent Variable 1",
      y = "Latent Variable 2",
      title = "Autoencoder"
    ) +
    scale_size(
      range = c(0.1, 15), 
      limits = c(1000, 1000000)
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
)  

# Saving cluster fasta
write_cluster_fasta(
  df = ae_features, 
  metagenome_fasta = fasta,
  path_out = paste0(p$misc$output_folder, "/bin_fasta_AE/")
)


