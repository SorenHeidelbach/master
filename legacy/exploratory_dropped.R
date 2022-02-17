
# Cluster preservation matrix
combinations <- expand_grid(x = unique(joined_features$type), y = unique(joined_features$cluster))
asd <- map2_dfr(combinations$x, combinations$y, function(x1, y1) {
  contigs <- joined_features %>% 
    filter(type == x1 & cluster == y1) %>% 
    pull(chr) 
  map2_dfc(combinations$x, combinations$y , function(x2, y2) {
    sum(contigs %in% (filter(joined_features, type == x2 & cluster == y2) %>% pull(chr)))
  })
}) %>% 
  setNames(paste0(combinations$x, " ", combinations$y)) %>% 
  mutate(from = paste0(combinations$x, " & ", combinations$y)) %>% 
  pivot_longer(cols = -from, names_to = "to") %>% 
  ggplot(aes(y = from, x = to, fill = value)) +
  geom_tile() + 
  scale_fill_steps(low = "white", 
                   high = "red",
                   limits = c(0, 60),
                   guide = guide_legend(title = "Shared\nContigs")) +
  geom_text(aes(label = ifelse(value == 0, "", value))) +
  geom_hline(yintercept = seq(0.5, 5*length(joined_features$cluster %>% unique) + 0.5, by = length(joined_features$cluster %>% unique))) +
  geom_vline(xintercept =  seq(0.5, 5*length(joined_features$cluster %>% unique) + 0.5, by = length(joined_features$cluster %>% unique))) +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 0),
        axis.text.y.left = element_text(hjust = 0)) +
  ggtitle("Cluster preservation with and without GC and coverage as features") +
  ylab("") +
  xlab("")

ggsave(
  filename = paste0(p$misc$output_folder, "AE_vs_tSNE_conservation_matrix.pdf"),  
  width = 8, height = 6, 
  asd
)





##-----------------------------------------------------------------------------
# PCA with coverage and GC
##-----------------------------------------------------------------------------

# With GC and coverage
pca_list_gc_cov <- MGM1_features %>% 
  left_join(coverage %>% select(chr, avg_cov.dataset_B, gc), by = c("contig" = "chr")) %>% 
  select(-contig) %>% 
  normalise_cols() %>% 
  pca_list()

pca_features_gc_cov <- MGM1_features %>% 
  left_join(coverage %>% select(chr, avg_cov.dataset_B, gc), by = c("contig" = "chr")) %>% 
  select(-contig) %>% 
  normalise_cols() %>% 
  #standardise_cols() %>% 
  pca() %>%
  mutate(
    cluster = hdbscan(cbind(PC1, PC2), minPts = 13)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = MGM1_features$contig) %>% 
  left_join(coverage)

ggsave(filename = "./figures/feature_exploration/PCA_ND-MGM1_hdbscan_with-gc-cov.pdf",  width = 8, height = 6,
       ggplot(pca_features_gc_cov, aes(x=PC1, y=PC2, fill=cluster, size = contig_length.x)) +
         theme_bw()+
         stat_density_2d(col = "gray60") +
         geom_point(aes(fill = cluster), shape = 21, alpha = 0.6) +
         geom_segment(data = pca_list_gc_cov$loading_vector, 
                      inherit.aes = FALSE,
                      aes(x = PC1*12, y = PC2*12),
                      xend = 0, yend = 0,
                      arrow = arrow(length = unit(0.01, "cm"))) +
         geom_text(data = filter(pca_list_gc_cov$loading_vector, (PC1^2 + PC2^2)^0.5 > 0.195),
                   inherit.aes = FALSE,
                   aes(label = label, x= PC1*12+0.5, y = PC2*12+0.5)
         ) +
         scale_size(range = c(1, 15), limits = c(25000, 900000)) +
         labs(x = paste0("PC1, explained variance: ", signif(100*pca_list_gc_cov$variance$var_explained[pca_list_gc_cov$variance$PC == "PC1"], 3), "%"), 
              y = paste0("PC2, explained variance: ", signif(100*pca_list_gc_cov$variance$var_explained[pca_list_gc_cov$variance$PC == "PC2"], 3), "%"), 
              title = "With Coverage and GC% features") +
         scale_size(guide = guide_legend(title = "Contig Length"))
)

##-----------------------------------------------------------------------------
# t-SNE with coverage and GC
##-----------------------------------------------------------------------------

# Feature transformation
tsne_features_gc_cov <- MGM1_features %>% 
  left_join(coverage %>% select(chr, avg_cov.dataset_B, gc), by = c("contig" = "chr")) %>% 
  select(-contig) %>% 
  normalise_cols() %>% 
  #standardise_cols() %>% 
  tsne() %>%
  mutate(
    cluster = hdbscan(cbind(V1, V2), minPts = 13)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = MGM1_features$contig) %>% 
  left_join(coverage)

# Plot
ggsave(filename = "./figures/feature_exploration/tsne_ND-MGM1_hdbscan_with-gc-cov.pdf",  width = 8, height = 6,
       ggplot(tsne_features_gc_cov, aes(x=V1, y=V2, fill=cluster, size = contig_length.x)) +
         stat_density_2d(col = "gray60") +
         geom_point(aes(fill = cluster), shape = 21, alpha = 0.6) +
         scale_size(range = c(1, 15), limits = c(25000, 900000)) +
         labs(x = "t-SNE dim 1", y = "t-SNE dim 2", title = "With Coverage and GC% features") +
         theme_bw()
)

# Cluster relationship matrix
joined_features <- rbind(
  mutate(tsne_features_gc_cov, gc_cov = "Meth only"),
  mutate(tsne_features, gc_cov = "M+GC+Cov")
)

combinations <- expand_grid(x = unique(joined_features$gc_cov), y = unique(joined_features$cluster))
asd <- map2_dfr(combinations$x, combinations$y, function(x1, y1) {
  contigs <- joined_features %>% 
    filter(gc_cov == x1 & cluster == y1) %>% 
    pull(chr) 
  map2_dfc(combinations$x, combinations$y , function(x2, y2) {
    sum(contigs %in% (filter(joined_features, gc_cov == x2 & cluster == y2) %>% pull(chr)))
  }
  )
}) %>% 
  setNames(paste0(combinations$x, " ", combinations$y)) %>% 
  mutate(from = paste0(combinations$x, " & ", combinations$y)) %>% 
  pivot_longer(cols = -from, names_to = "to") %>% 
  ggplot(aes(y = from, x = to, fill = value)) +
  geom_tile() + 
  scale_fill_steps(low = "white", 
                   high = "red",
                   limits = c(0, 60),
                   guide = guide_legend(title = "Shared\nContigs")) +
  geom_text(aes(label = ifelse(value == 0, "", value))) +
  geom_hline(yintercept = seq(0.5, 5*length(joined_features$cluster %>% unique) + 0.5, by = length(joined_features$cluster %>% unique))) +
  geom_vline(xintercept =  seq(0.5, 5*length(joined_features$cluster %>% unique) + 0.5, by = length(joined_features$cluster %>% unique))) +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 0),
        axis.text.y.left = element_text(hjust = 0)) +
  ggtitle("Cluster preservation with and without GC and coverage as features") +
  ylab("")+
  xlab(" ")
ggsave(filename = glue("./figures/feature_exploration/MGM1_cluster_conservation_matrix.pdf"), plot = asd)



##-----------------------------------------------------------------------------
# PCoA plotting
#c("euclidean", "minkowski", "maximum", "manhattan", "canberra", "binary") %>% 
#  map(function(method) {
#    MGM1_features %>% 
#      left_join(coverage %>% select(chr, avg_cov.dataset_B, gc), by = c("contig" = "chr")) %>% 
#      mutate(
#        avg_cov.dataset_B = log(avg_cov.dataset_B)
#      )
#      select(-contig) %>% 
#      standardise_cols() %>% 
#      dist(method = method) %>% 
#      ape::pcoa() %>% 
#      `[[`("vectors") %>% 
#      as.data.frame() %>% 
#      mutate(
#        cluster = hdbscan(cbind(Axis.1, Axis.2), minPts = 13)$cluster,
#        cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
#      ) %>%
#      cbind(chr = MGM1_features$contig) %>% 
#      left_join(coverage) %>% 
#      ggplot(aes(x = Axis.1, y = Axis.2, size = contig_length.x, fill = cluster)) +
#      theme_bw() +
#      geom_point(shape = 21, color = "black", alpha = 0.5) +
#      ggtitle(paste0("Distance method: ", method)) +
#      labs(subtitle = "With GC% and coverage")
#)









# Checking position correlation within different motifs
MGM1_motif %>% 
  unnest(cols = c(data)) %>% 
  filter(motif %in% c("GATC", "GGGCTC") & contig %in% pull(filter(tsne_features, cluster == "3"), chr)) %>% 
  mutate(feature = paste0(motif, "_",distance_motif)) %>% 
  select(contig, feature, dist_score) %>% 
  pivot_wider(values_from = dist_score, names_from = feature) %>% 
  select(-contig) %>% 
  mutate_all(
    function(x) ifelse(is.na(x), 0, x)
  ) %>% 
  normalise_cols() %>% 
  mutate_all(function(x) log(x+1)) %>% 
  chart.Correlation(histogram=TRUE, pch=19, method = "spearman")


