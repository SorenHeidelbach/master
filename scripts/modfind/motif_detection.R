#!/usr/bin/Rscript
.libPaths(c("/shared-nfs/SH/code/Rlib", .libPaths()))
setwd(here::here())
require("argparse")

########################################################################################################
# Argument parser
########################################################################################################

# Read in arguments
parser <- ArgumentParser()
parser$add_argument("-s", "--style", nargs = 1, help = "Style embedding of Dual-AAE")
parser$add_argument("-c", "--cat", nargs = 1, help = "Categerical embedding of Dual-AAE")
parser$add_argument("-a", "--assembly", nargs = 1, help = "Assembly with contigs")
parser$add_argument("-o", "--out", nargs = 1, help = "Output folder")
parser$add_argument("-f", "--features", nargs = 1, help = "Input features to Dual-AAE")
parser$add_argument("-m", "--metainfo", nargs = 1, help = "Metainfo outputted from preprocess")
parser$add_argument("--metabat", nargs = 1, help = "tsv file with one column of contigs and another with bins from metabat")
parser$add_argument("--assembly_info", nargs = 1, help = "Assembly_info file from assembly")
parser$add_argument("--GTDBTk", nargs = 1, help = "Summary file from GTDBTk")
parser$add_argument("--run_umap", action="store_true", help = "Include UMAP embedding in analysis (Can be compute intensive)")

parser$add_argument("--train_loss", nargs = 1, help = "Metainfo outputted from preprocess")
parser$add_argument("--validation_loss", nargs = 1, help = "Metainfo outputted from preprocess")

parser$add_argument("--min_cluster_size", default = 50, help = "Minimum number of event in cluster to included")
parser$add_argument("-t", "--threads", default = 10, type = "integer", help = "Minimum number of event in cluster to included")

args <- parser$parse_args()


########################################################################################################
# Argument
########################################################################################################

# Check required arguments are present
if (is.null(args$style)) {
  stop(sprintf("Required argument \"-s, --style\" is missing"))
}
if (is.null(args$cat)) {
  stop(sprintf("Required argument(\"-c, --cat\") is missing"))
}
if (is.null(args$assembly)) {
  stop(sprintf("Required argument(\"-a, --assembly\") is missing"))
}

dir.create(paste0(args$out, "/DAAE/plots"), recursive = T, showWarnings = FALSE)
dir.create(paste0(args$out, "/UMAP/plots"), recursive = T, showWarnings = FALSE)


########################################################################################################
# Load packages
########################################################################################################

pacman::p_load("data.table", "seqinr", "msa", "dbscan", "Rtsne", "dplyr",
               "ggplot2", "tidyr", "jsonlite", "stringr", "umap", "ggforce", "mice", "visdat")

source("/shared-nfs/SH/code/scripts/functions.R")


# Save setting for detection
write_json(args, paste0(args$out, "/settings_motif_detection.json"), pretty = TRUE)


########################################################################################################
# Debugging parameters
########################################################################################################
#
# args = list()
# args$style <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/DAAE/embedding_style_epoch_60.tsv'
# args$categorical <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/DAAE/embedding_categorical_epoch_60.tsv'
# args$out <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6'
# args$features <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/event_features.tsv'
# args$metainfo <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/metainfo.tsv'
# args$threads <- 50
# args$assembly <- '/shared-nfs/SH/samples/zymo/assembly_100x.fasta'
# args$min_cluster_size <- 50
# args$validation_loss <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/DAAE/loss_validation_epoch_60.tsv'
# args$train_loss <- '/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_6/DAAE/loss_training_epoch_60.tsv'
# args$metabat <- '/shared-nfs/SH/samples/zymo/metabat2/bin_membership.tsv'
# args$assembly_info <- '/shared-nfs/SH/samples/zymo/NAT/assembly_filt_100x_trim/flye/assembly_info.txt'
# args$GTDBTk <- '/shared-nfs/SH/samples/zymo/metabat2/GTDBTk/classify/gtdbtk.bac120.summary.tsv'

########################################################################################################
# Loss plot
########################################################################################################

if (!(is.null(args$train_loss) & is.null(args$validation_loss))) {
  loss_val <- fread(args$train_loss) %>%
    cbind(type = "train", epoch = 1:nrow(.))
  loss_train <- fread(args$validation_loss) %>%
    cbind(type = "validation", epoch = 1:nrow(.))
  loss <- rbind(loss_val, loss_train) %>%
    melt(id.vars = c("type", "epoch"))
  ggsave(
    filename = paste0(args$out, "/DAAE/plots/loss.png"),
    device = "png",
    height = unit(7, "cm"),
    width = unit(10, "cm"),
    ggplot(loss, aes(x = epoch, y = value, col = type)) +
      geom_line() + facet_wrap(~variable, scales = "free") +
      theme_bw() +
      labs(x = "Epoch", y = "Loss", title = "Losses during training of Dual-AAE") +
      guides(col = guide_legend(title = NULL)))
}


########################################################################################################
# Load data
########################################################################################################

print("Reading in data")
event <- fread(args$style) %>%
  setnames(paste0("style_", 1:ncol(.)))
event <- event %>%
  cbind(
    fread(args$metainfo)[1:nrow(.)]
  )
event <- event %>%
  cbind(
    fread(args$features)[1:nrow(.)]
  )
event <- event %>%
  cbind(
    fread(args$cat) %>%
      setnames(paste0("group_", 1:ncol(.))) %>%
      apply(MARGIN = 1, FUN = function(x) {
        i = which(x == max(x))[1]
        c(names(x)[i], x[i])
      }) %>%
      t() %>%
      data.frame(stringsAsFactors = FALSE) %>%
      `names<-`(c("DAAE_group", "DAAE_group_prob")) %>%
      mutate(
        DAAE_group_prob = as.numeric(DAAE_group_prob),
        DAAE_group = ifelse(DAAE_group_prob < 0.05, NA, DAAE_group)
        )
  )
# Getting sequence information
assembly <- read.fasta(args$assembly)

event <- event %>%
  left_join(
    lapply(assembly, length) %>% unlist %>% data.frame(length = ., contig = names(.)),
    by = "contig"
  )

########################################################################################################
# UMAP on features
########################################################################################################
if (args$run_umap) {
  print("UMAP on all features")
  if (file.access(paste0(args$out, "/UMAP/embedding.RDS"), 4) == -1) {
    umap_all <- umap(select(event, starts_with(c("rel", "mean"))), verbose = TRUE, n_epochs = 150)
    saveRDS(umap_all, file = paste0(args$out, "/UMAP/embedding.RDS"))
  } else {
    umap_all <- readRDS(paste0(args$out, "/UMAP/embedding.RDS"))
  }

  event <- event %>%
    cbind(
      umap_all$layout %>%
        data.frame() %>%
        setnames(c("UMAP1", "UMAP2"))
    )

  # DBSCAN
  print("DBSCAN on UMAP embedding")
  dbscan_umap <- dbscan(select(event, starts_with("UMAP")), minPts = 50, eps = 0.05)

  # Add clusters to df
  event <- event %>%
    cbind(UMAP_DBSCAN = dbscan_umap$cluster)

  # Find biggest clusters for later viz
  event <- event %>%
    mutate(
      top_n_dbscan = UMAP_DBSCAN %in% (
        event %>%
          group_by(UMAP_DBSCAN) %>%
          summarise(count = n()) %>%
          filter(count %in% tail(sort((count)), 20)) %>%
          pull(UMAP_DBSCAN))
    )
}

########################################################################################################
# Dual-AAE
########################################################################################################

# Dual-AAE with clustering
print("Dual-AAE plots")

##
# Style on xy and colored by categorical
##
print("  1")
print(event)
ggsave(
  filename = paste0(args$out, "/DAAE/plots/style_with_categorical.png"),
  device = "png",
  width = unit(10, "cm"),
  height = unit(7, "cm"),
  event %>%
    group_by(DAAE_group) %>%
    mutate(DAAE_group = ifelse(n() > args$min_cluster_size, DAAE_group, NA)) %>%
    ggplot(aes(x = style_1, y = style_2, col = DAAE_group)) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw() +
    stat_density_2d(geom = "polygon",
                    contour = TRUE,
                    aes(fill = after_stat(level)),
                    colour = "black",
                    alpha = 0.3,
                    bins = 5) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    labs(
      x = "Style dimension 1",
      y = "Style dimension 2",
      title = "Dual-AAE clustering") +
    guides(
      fill = guide_legend(title = "Event \nDensity"),
      col = guide_legend(title = "Categorical"))
  )
##
# Number of events in each group
##
print("  2")
ggsave(
  paste0(args$out, "/DAAE/plots/n_events_in_categoricals.png"),
  device = "png",
    event %>%
    group_by(DAAE_group) %>%
    summarise(count = n()) %>%
    mutate(DAAE_group = factor(DAAE_group, levels = unique(DAAE_group[order(count, decreasing = TRUE)]))) %>%
    ggplot(aes(x = DAAE_group, y = count)) +
    geom_bar(stat = "identity") +
    scale_y_log10() +
    theme(axis.text.x.bottom = element_text(angle = 90)) +
    geom_hline(yintercept = args$min_cluster_size))

########################################################################################################
# Dual-AAe, with assembly info and binning
########################################################################################################
if (!(is.null(args$GTDBTk) & is.null(args$metabat))) {
  # Evaluate in context of Metabat
  metabat <- fread(args$metabat)
  GTDBTk <- fread(args$GTDBTk)

  ggsave(
    paste0(args$out, "/DAAE/plots/groups_metabat_species.png"),
    width = unit(12, "cm"),
    height = unit(8, "cm"),
    device = "png",
    event %>%
      group_by(contig, DAAE_group) %>%
      summarise(count = n()) %>%
      left_join(metabat, by = "contig") %>%
      left_join(GTDBTk, by = c(bin = "user_genome")) %>%
      mutate(species = str_extract(classification, "s__.*")) %>%
      ungroup() %>%
      group_by(DAAE_group) %>%
      filter(max(count) > 100) %>%
      ungroup() %>%
      mutate(DAAE_group = factor(DAAE_group, levels = unique(DAAE_group[order(count)]))) %>%
      ggplot(aes(x = DAAE_group, y = count)) + geom_bar(stat = "identity") + facet_wrap(~species,
      scales = "free") + theme(axis.text.x.bottom = element_text(angle = 90, size = 4)))
}

if (args$run_umap) {
  print("UMAP plots")
  # UMAP plot
  ggsave(
    filename = paste0(args$out, "/UMAP/plots/scatter.png"),
    device = "png",
    height = unit(6, "cm"),
    width = unit(9, "cm"),
    ggplot(event, aes(x = UMAP1, y = UMAP2)) +
      geom_point(alpha = 0.01, size = 0.6) +
      stat_density_2d(geom = "polygon", contour = TRUE,
                      aes(fill = after_stat(level)), colour = "black",
                      bins = 5) +
      scale_fill_distiller(palette = "Reds", direction = 1) +
      theme_bw() +
      labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP embedding of event features") +
      guides(fill = guide_legend(title = "Event\nDensity"))
  )

  # UMAP with DBSCAN
  ggsave(
    filename = paste0(args$out, "/UMAP/plots/scatter_dbscan.png"),
    device = "png",
    height = unit(6, "cm"),
    width = unit(13, "cm"),
    ggplot(event, aes(x = UMAP1, y = UMAP2, fill = as.character(UMAP_DBSCAN))) +
      geom_point(alpha = 0.01, size = 1.5, shape = 21) +
      theme_bw() +
      geom_text(
        data = filter(event, top_n_dbscan) %>%
          group_by(UMAP_DBSCAN) %>%
          summarise(
            x_mean = mean(UMAP1),
            y_mean = mean(UMAP2)
            ),
        aes(x = x_mean, y = y_mean, label = UMAP_DBSCAN)) +
      labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP embedding of event features with DBSCAN clusters") +
      guides(fill = guide_legend(title = "Cluster", override.aes = list(alpha=1))))

  # UMAP with DAAE categoricals
  ggsave(
    filename = paste0(args$out, "/UMAP/plots/scatter_DAAE.png"),
    device = "png",
    height = unit(6, "cm"),
    width = unit(13, "cm"),
    event %>%
    group_by(DAAE_group) %>%
    mutate(DAAE_group = ifelse(n() > args$min_cluster_size, DAAE_group, NA)) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, fill = DAAE_group)) +
    geom_point(alpha = 0.1, size = 1.5, shape = 21) +
    theme_bw() +
    labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP embedding of event features with DAAE clusters") +
    guides(fill = guide_legend(title = "Cluster", override.aes = list(alpha=1))))
}
########################################################################################################
# Cluster consensus sequence
########################################################################################################

MSA <- function(cluster_col){
  #
  #   Get sequence of event frames and perform multiple sequence alignment of
  #   each cluster.
  #   Return list with sequences, MSA objects and data frame
  #

  out <- list()
  out$seq <- mapply(function(contig, pos, DAAE_group) {
    x = paste0(assembly[[contig]][(max(1, pos - 3)):(pos + 3)], collapse = "")
    names(x) = DAAE_group
    return(x)
  }, event$contig, event$position, event[[cluster_col]])

  # Remove contig from group name
  names(out$seq) <- gsub(".*\\.", "", names(out$seq))

  # Multiple sequence alignment of each group
  out$seq2 <- lapply(unique(names(out$seq)), function(DAAE_group) {
    seq <- out$seq[names(out$seq) %in% DAAE_group] %>%
      DNAStringSet()
  })
  out$seq2 <- out$seq2[lapply(out$seq2, length) > 1]
  out$MSA <- out$seq2 %>%
    lapply(function(x) {
      if (length(x) > 500) {
        x = x[sample.int(n = length(x), size = 500)]
      }
      if (length(x) > 1) {
        x %>%
          msa("ClustalW", type = "dna")
      }
    })
  out$MSA_consensus <- out$MSA %>%
    lapply(msaConsensusSequence, type = "upperlower", thresh = c(100,90))
  out$df <- data.frame(
    cluster = out$MSA %>%
      lapply(function(x) x@unmasked@ranges@NAMES[1]) %>%
      unlist(),
    sequence = out$MSA_consensus %>%  unlist() %>%
      str_remove_all("-"),
    n_seqs = lapply(out$seq2, length) %>% unlist()
    )

  out
}
MSA_DAAE <- MSA("DAAE_group")
saveRDS(MSA_DAAE, paste0(args$out, "/DAAE/msa.RDS"))
if (args$run_umap) {
  MSA_UMAP <- MSA("UMAP_DBSCAN")
  saveRDS(MSA_UMAP, paste0(args$out, "/UMAP/msa.RDS"))
}

########################################################################################################
# Contig binning using new features (UMAP embeddgin of features)
########################################################################################################

contig_events <- function(event_df, group_col, minPts = 15){
  #
  #   Convert event clusters from embedding+clustering into features for contig
  #   binning. Use the column name in the df to specify which modification groups
  #   to use for binning
  #

  # Calculate relative numer of event clusters for each contig
  binning_features <- event_df %>%
    group_by(contig, !!sym(group_col), length) %>%
    summarise(count = n()) %>%
    group_by(contig) %>%
    mutate(
      rel_count = count/sum(count),
      cluster = paste0("cluster_", !!sym(group_col))) %>%
    ungroup() %>%
    select(contig, length, cluster, rel_count) %>%
    pivot_wider(names_from = cluster, values_from = rel_count)
  missing <- md.pattern(binning_features, rotate.names = TRUE)

  # Assign missing identifier to NAs
  binning_features_impute <- binning_features
  binning_features_impute[is.na(binning_features_impute)] <- -1
  # UMAP embedding of relative nubmer event clusters
  contig_umap <- binning_features_impute %>%
    select(-contig, -length) %>%
    umap()

  # HDBSCAN in embedding for contig clustering
  binned_contigs <- contig_umap %>%
    `[[`(1) %>%
    hdbscan(minPts = minPts) %>%
    `[[`("cluster") %>%
    data.frame(
      cluster = .,
      contig_umap$layout,
      binning_features
    )
  return(
    list(
      "UMAP" = contig_umap,
      "binned_contigs" = binned_contigs,
      "data_missingness" = missing
    )
  )
}


plot_contig_events <- function(clustered_contigs_df){
  #
  #   Scatter plot of UMAP embedding with HDBSCAN clusters of modification
  #   features
  #

  df <- clustered_contigs_df
  df %>%
    left_join(
      event %>%
        group_by(contig) %>%
        filter(DAAE_group != "NA") %>%
        summarise(count = n())
    ) %>%
    ggplot(aes(x = X1, y = X2, fill = as.character(cluster))) +
    geom_point(aes(size = length, alpha = count), shape = 21) +
    geom_label(
      data = (df %>%
                group_by(cluster) %>%
                summarise(
                  size = signif(sum(length)/1e6, 3),
                  ymean = mean(X2),
                  xmean = mean(X1)) %>%
                filter(size > 1) %>%
                mutate(size = paste0(size, " Mb"))),
      mapping = aes(label = size, x = xmean, y = ymean),
      fill = "gray99",
      label.padding = unit(0.15, "lines")
    ) +
    labs(x = "UMAP1", y = "UMAP2", title = "UMAP embedding of relative event cluster frequency") +
    scale_size(trans = "log10") +
    scale_alpha(trans = "log10") +
    theme_bw() +
    theme(
      text = element_text(size = 18),
      axis.text = element_text(size = 10)
    ) +
    guides(
      fill = guide_legend(title = "Bin"),
      size = guide_legend(title = "Contig \nLength"),
      alpha = guide_legend(title = "Number of \nEvents")
    )
}

plot_categorical <- function(clustered_contigs_df, cluster_col) {
  #
  #   Number of events in modification clusters, facetted by bins
  #

  event %>%
    left_join(
      clustered_contigs_df %>%
        select(contig, cluster)
    ) %>%
    filter(!is.na(!!sym(cluster_col))) %>%
    group_by(cluster, !!sym(cluster_col)) %>%
    summarise(count = n()) %>%
    mutate(cluster = paste0("Bin ", cluster)) %>%
    ggplot(aes(x = as.character(!!sym(cluster_col)), y = count)) +
    geom_bar(stat = "identity") +
    facet_wrap(~cluster, scales = "free_y") +
    #scale_y_log10() +
    theme_minimal() +
    theme(
      text = element_text(size = 18),
      axis.text = element_text(angle = 90, hjust = 1, size = 10)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Modification Event Groups", y = "Number of Events")
}

plot_missing <- function(missing_df){
  missing_df %>%
    select(starts_with("cluster_")) %>%
    vis_miss(warn_large_data = FALSE) +
    labs(y = "Contigs") +
    theme(
      text = element_text(size = 18),
      axis.text = element_text(size = 10),
      axis.text.x.top = element_text(angle = 90)
    )
}

contig_binning <- function(out_dir, cluster_col){
  binning <- contig_events(event, cluster_col)
  ggsave(
    paste0(args$out, "/", out_dir, "/plots/contig_binning_", cluster_col,".png"),
    width = unit(12, "cm"),
    height = unit(8, "cm"),
    device = "png",
    plot_contig_events(binning$binned_contigs)
  )
  ggsave(
    paste0(args$out, "/", out_dir, "/plots/event_group_count_", cluster_col,".png"),
    width = unit(12, "cm"),
    height = unit(8, "cm"),
    device = "png",
    plot_categorical(binning$binned_contigs, cluster_col)
  )
  ggsave(
    paste0(args$out, "/", out_dir, "/plots/feature_missingness_", cluster_col,".png"),
    width = unit(12, "cm"),
    height = unit(8, "cm"),
    device = "png",
    plot_missing(binning$binned_contigs)
  )
  return(binning)
}

# DAAE
DAAE_binned <- contig_binning("DAAE", "DAAE_group")
write_cluster_fasta(
  df = DAAE_binned$binned_contigs,
  metagenome_fasta = assembly,
  path_out = paste0(args$out, "/DAAE/bins/")
)

# UMAP
if (args$run_umap) {
  UMAP_binned <- contig_binning("UMAP", "UMAP_DBSCAN")
  write_cluster_fasta(
    df = UMAP_binned$binned_contigs,
    metagenome_fasta = assembly,
    path_out = paste0(args$out, "/UMAP/bins/")
  )
}
