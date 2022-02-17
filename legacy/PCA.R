
pca_list <- features %>% 
  select(-contig) %>% 
  normalise_cols() %>% 
  pca_list()

pca_features <- features %>% 
  select(-contig) %>% 
  normalise_cols() %>% 
  pca() %>%
  mutate(
    cluster = hdbscan(cbind(PC1, PC2), minPts = p$hdbscan$minPts)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = features$contig) %>% 
  left_join(coverage)
write_csv(pca_features, file = paste0(p$misc$output_folder, "/features_PCA.csv"))
ggsave(
  filename = paste0(p$misc$output_folder, "/cluster_PCA.pdf"), 
  width = 8, height = 6,
  ggplot(pca_features, aes(x=PC1, y=PC2, fill=cluster)) +
    theme_bw() +
    geom_point(
      aes(size = contig_length.x), 
      alpha = 0.4, 
      shape = 21
    ) +
    geom_segment(
      data = pca_list$loading_vector, 
      inherit.aes = FALSE,
      aes(x = PC1*12, y = PC2*12),
      xend = 0, yend = 0,
      arrow = arrow(length = unit(0.01, "cm"))
    ) +
    geom_text(
      data = filter(pca_list$loading_vector, (PC1^2 + PC2^2)^0.5 > 0.195),
      inherit.aes = FALSE,
      aes(label = label, x= PC1*12+0.5, y = PC2*12+0.5)
    ) +
    scale_size(
      range = c(1, 15), 
      limits = c(25000, 900000)
    ) +
    labs(
      x = paste0("PC1, explained variance: ", signif(100*pca_list$variance$var_explained[pca_list$variance$PC == "PC1"], 3), "%"), 
      y = paste0("PC2, explained variance: ", signif(100*pca_list$variance$var_explained[pca_list$variance$PC == "PC2"], 3), "%")
    ) +
    scale_size(guide = guide_legend(title = "Contig Length"))
)

ggsave(
  filename = paste0(p$misc$output_folder, "/PCA_scree.pdf"),  
  width = 8, height = 6,
  pca_list$variance %>% 
    filter(var_explained > 0.005) %>% 
    ggplot(aes(x = PC, y = var_explained)) +
    geom_bar(stat = "identity") +
    ylab("Variance Explained") +
    theme_bw() +
    theme(axis.text.x.bottom = element_text(angle = 90, hjust = 0)) 
)
# Saving cluster fasta
write_cluster_fasta(
  df = pca_features, 
  metagenome_fasta = fasta,
  path_out = paste0(p$misc$output_folder, "/bin_fasta_PCA/")
)