
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


###########################################################################
## Nanodisco feature engineering
###############################################################################

##-----------------------------------------------------------------------------
# Loading NanoDisco profile
load_ND_profile <- function(path,
                            min_contig_size = 25000) {
  # Load profile
  message("Loading profile")
  profile <- readRDS(path)
  profile <- as.data.table(profile)
  message("  Metainformation")
  # Contig meta information
  contigs <- profile %>%
    attributes() %>%
    `[[`("contig_coverage") %>%
    filter(contig_length > min_contig_size)  %>%
    arrange(-contig_length) %>%
    mutate(size_order = as.numeric(row.names(.))) %>%
    mutate(contig=chr) %>%
    select(-chr)
  
  message("  Removing small contigs")
  # Combine contig information with profile data frame
  profile[
    contigs, on = "contig", (names(contigs)[-7]) := mget(paste0("i.", names(contigs)[-7]))
    ][
      contig_length >= min_contig_size
      ]
  return(profile)
}

##-----------------------------------------------------------------------------
# Profile into features
get_ND_profile_features <- function(profile,
                                    min_nb_occurrence = 10,
                                    min_dist_score = 1.5,
                                    min_n = 10,
                                    filter_methods = "max",
                                    missing_value = 0,
                                    select_features = "all"){
  # Loading nested profile dataframe
  profile
  message("Filtering")
  if(select_features != "all"){
    profile %>%
      filter(paste0(motif, "_", distance_motif) %in% select_features)
  }
  
  # Removing feature vectors that doesn't contain any informative values
  if(filter_methods == "max"){
    message("  Keeping all values in feature vector if just one of the values are above minimum")
    ND_features <- profile %>%
      group_by(motif, distance_motif) %>%
      filter(max(nb_occurrence) > min_nb_occurrence & max(dist_score) > min_dist_score)
  } else {
    message("  Choose valid filtering method")
    stop()
  }
  
  message("Indentifying missing feature values")
  # pviot longer to get each feature in a column
  ND_features <- ND_features %>%
    ungroup() %>%
    mutate(
      feature_name = paste0(motif, "_", distance_motif)
    ) %>%
    select(contig, feature_name, dist_score) %>%
    pivot_wider(values_from = dist_score, names_from = feature_name)
  
  ND_features <- ND_features %>% pivot_longer(!contig, names_to = "feature_name", values_to = "dist_score")
  
  missingness <- sum(is.na(ND_features))/(nrow(ND_features)*ncol(ND_features))
  message(paste0("  Percent of missing values: ", missingness))
  # Setting missing values according to specified method
  contigs <- ND_features$contig
  if (missing_value %in% "mean") {
    message("  Mising values are filled with feature vector mean")
    ND_features <- ND_features %>%
      group_by(feature_name) %>%
      mutate(
        dist_score = ifelse(is.na(dist_score), mean(dist_score, na.rm = TRUE), dist_score)
      )
  } else if (is.numeric(missing_value)) {
    message(paste0("  Missing values are replaced with ", missing_value))
    ND_features <- ND_features %>%
      group_by(feature_name) %>%
      mutate(
        dist_score = ifelse(is.na(dist_score), missing_value, dist_score)
      )
  } else if (missing_value == "NA") {
    message("  Missing values are kept as NA")
    # returns dataframe with NA values
  } else {
    message("  Please set a valid 'missing_value' ('mean' or replacement number or 'NA')")
    stop()
  }
  ND_features <- ND_features %>%
    pivot_wider(values_from = dist_score, names_from = feature_name)
  message("  Missing values treated")
  return(ND_features)
}



###############################################################################
## Ordination/dimontionality reduction
###############################################################################
# All functions take in a feature matrix with an ID column and a metadata dataframe

##-----------------------------------------------------------------------------
# PCA

pca <- function(features, metadata = NA, nPC = 2) {
  df <- prcomp(features, center = TRUE, scale. = TRUE) %>%
    `[[`(5) %>%
    as.data.frame() %>%
    select(paste0("PC", 1:nPC)) %>%
    mutate(contig = features$contig)
  if (!is.na(metadata)) {
    attr(df, "metadata") <- metadata
  }
  df
}

# PCA relevant information
pca_list <- function(df){
  # PCA object
  pca_object <- df %>%
    prcomp()
  
  # PCA explained variance
  pca_explained_var <- data.frame(
    PC = (pca_object$x %>% colnames()),
    var_explained=(pca_object$sdev)^2/sum((pca_object$sdev)^2)
  ) %>%
    mutate(PC = factor(PC, levels = PC[order(var_explained, decreasing = TRUE)])) %>%
    arrange(-var_explained)
  
  # PCA loading vector for each PC
  pca_loading_vector <- pca_object$rotation %>%
    data.frame() %>%
    mutate(label = rownames(.))
  return(list(object = pca_object, variance = pca_explained_var, loading_vector = pca_loading_vector))
}

##-----------------------------------------------------------------------------
# t-SNE

tsne <- function(features, metadata = NA, seed = 1, ...){
  set.seed(seed)
  tsne_obj <- Rtsne(as.matrix(features), ...)
  df <- tsne_obj$Y %>%
    as.data.frame()
  if (!is.na(metadata)) {
    attr(df, "metadata") <- metadata
  }
  df
}
tsne_obj <- function(features,
                     theta = 0.5,
                     perplexity = 5,
                     max_iter = 1000,
                     metadata = NA,
                     seed = 1){
  set.seed(seed)
  tsne_obj <- Rtsne(as.matrix(features),
                    theta = theta,
                    perplexity = perplexity,
                    max_iter = max_iter)
  df <- tsne_obj$Y %>%
    as.data.frame()
  if (!is.na(metadata)) {
    attr(df, "metadata") <- metadata
  }
  tsne_obj
}

##-----------------------------------------------------------------------------
# autoencoder

autoencoder <- function(features,
                        bottleneck_nodes = 2,
                        intermediate_nodes = 6,
                        epochs = 1000,
                        activator = "relu"){
  message("Defining model")
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = intermediate_nodes,
                activation = activator,
                input_shape = ncol(features)) %>%
    layer_batch_normalization() %>%
    layer_dense(units = intermediate_nodes,
                activation = activator) %>%
    layer_batch_normalization() %>%
    layer_dense(units = bottleneck_nodes,
                activation = activator,
                name = "bottleneck") %>%
    layer_batch_normalization() %>%
    layer_dense(units = intermediate_nodes,
                activation = activator) %>%
    layer_batch_normalization() %>%
    layer_dense(units = intermediate_nodes,
                activation = activator) %>%
    layer_batch_normalization() %>%
    layer_dense(units = ncol(features))
  
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = "adam"
  )
  
  message("Training model")
  features_train <- as.matrix(features)
  model %>% fit(
    x = features_train,
    y = features_train,
    epochs = epochs,
    verbose = 1
  )
  
  message("Extracting latent variables")
  # extract the bottleneck layer
  bottleneck_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
  
  latent_features <- predict(bottleneck_layer_model, features_train) %>%
    as.data.frame()
}


###############################################################################
##
###############################################################################

cluster_plot <- function(ggplot_obj, savepath,
                         xlab = "Variable 1", ylab = "Variable 2",
                         dim = c(8, 6)){
  ggplot_obj <- ggplot_obj + geom_density_2d(contour = TRUE, col = "gray60") +
    geom_point(shape = 21, alpha = 0.7) +
    theme_bw() +
    scale_size(
      range = c(0.1, 15),
      limits = c(1000, 1000000)
    ) +
    guides(
      fill = guide_legend(title="Cluster"),
      size = guide_legend(title = "Contig Size", override.aes = list(col=NA, fill = "black"))
    ) +
    
    labs(
      x = xlab,
      y = ylab
    )
  ggsave(ggplot_obj, filename = savepath, width = dim[1], height = dim[2])
}

