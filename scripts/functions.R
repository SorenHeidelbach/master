
# Setting costum library path
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
if (!require("pacman")) install.packages("pacman", lib = "/shared-nfs/SH/code/Rlib")
if (!require("BiocManager")) install.packages("BiocManager", lib = "/shared-nfs/SH/code/Rlib")

# Loading libraries with pacman
pacman::p_load(
  "ggplot2",
  "tidyr",
  "dplyr",
  "Rtsne",
  "purrr",
  "dbscan",
  "glue",
  "seqinr",
  "jsonlite"
)
##-----------------------------------------------------------------------------
# Calculate GC% from fasta file path
get_GC_percent <-  function(fasta, min_size = 0) {
  fasta %>%
    names() %>%
    map_df(
      function(contig) {
        seq <- fasta %>%
          `[[`(contig) %>%
          c() %in% c("g", "c", "G", "C")
        size = length(seq)
        gc <- sum(seq)/size
        data.frame(contig = contig, gc = gc, contig_length = size)
      }
    ) %>%
    filter(contig_length >= min_size)
}
##-----------------------------------------------------------------------------
# Normalise values in columns of data frame
normalise_cols <- function(numeric_dataframe){
  numeric_dataframe %>%
    apply(MARGIN = 2, function(x) (x-min(x)) / (max(x) - min(x))) %>%
    as.data.frame()
}
##-----------------------------------------------------------------------------
# Standardise values in columns of data frame
standardise_cols <- function(numeric_dataframe){
  numeric_dataframe %>%
    apply(MARGIN = 2, function(x) (x - mean(x))/sd(x)) %>%
    as.data.frame()
}

###############################################################################
## Misc functions
###############################################################################
# Read GFA file from assembly
read_gfa <- function(gfa_path){
  system(paste0("awk '/^P/ {print $2,$3}' ", gfa_path), intern = TRUE ) %>%
    data.frame() %>%
    separate(col = ".", into = c("contig", "edge"), sep = " ") %>%
    mutate(
      edge = str_remove(edge, "\\+")
    ) %>%
    left_join(
      system(paste0("awk '/^S/ {print $2,$4}' ", gfa_path), intern = TRUE ) %>%
        data.frame() %>%
        separate(col = ".", into = c("edge", "cov"), sep = " ") %>%
        mutate(cov = str_remove(cov, "dp\\:i\\:"))
    )
}
# Saving cluster
write_cluster_fasta <-  function(df,
                    metagenome_fasta,
                    path_out) {
  dir.create(path_out)
  df$cluster %>%
    unique()  %>%
    map(function(cluster_val){
      cluster_contigs <- df %>%
        filter(cluster == cluster_val) %>%
        pull(contig)
      seq <- metagenome_fasta[names(metagenome_fasta) %in% cluster_contigs]
      write.fasta(seq, names = cluster_contigs ,file.out = paste0(path_out, "cluster_", cluster_val, ".fasta"))
    })
  return("Cluster sequences written")
}

# Read in wierd checkm tsv files
read_checkm <- function(path){
  bin_stats_raw <- read_delim(path, col_names = FALSE )
  bin_header <- bin_stats_raw %>%
    separate(X2, into = as.character(1:(str_count(.$X2, "\\:"))), "[,]") %>%
    `[`(1,) %>%
    mutate(across(c(-X1),
                  function(x) {
                    y <- str_extract(x, "(?![\\{\\#]).*(?=\\:)")
                    y <- str_remove_all(y, "[\\'\\#]")
                    y <- str_trim(y)
                    y <- str_replace_all(y, " ", "_")
                  })) %>% c() %>%
    unlist %>%  c()
  bin_header[1] <- "bin"

  bin_stats <- bin_stats_raw %>%
    separate(X2, into = as.character(1:(str_count(.$X2, "\\:"))), "[,]") %>%
    mutate(across(c(-X1),
                  function(x) str_extract(x, "(?<=\\:)[^\\}]*"))) %>%
    set_names(bin_header)
}

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

tsne <- function(features,
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
