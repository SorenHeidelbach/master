#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
source("./code/scripts/functions.R")
pacman::p_load(
  "data.table",
  "purrr",
  "stringr"
)

# Read in parameters
p <- jsonlite::fromJSON("./code/params.json")

# Loading motifs
motif <- load_ND_profile(
  p$misc$ND_profile_path,
  min_contig_size = p$feature$min_contig_size
)


# Converting motif information to features matrix
features <- get_ND_profile_features(
  motif,
  filter_methods = p$feature$filter_method,
  min_dist_score = p$feature$min_dist_score,
  min_nb_occurrence = p$feature$min_nb_occurrence,
  missing_value = p$feature$missing_value
)
features_norm <- features %>%
  `[`(,2:ncol(.))

features_norm <- (features_norm - min(features_norm))/max(features_norm)
features_norm <- (features_norm - min(features_norm))/max(features_norm)
features_norm <- -1 + 2*features_norm
# Output folder
dir.create(p$misc$output_folder)

# Saving feature dataframe
write_csv(features, file = paste0(p$misc$output_folder, "/features.csv"))
write_csv(features_norm, file = paste0(p$misc$output_folder, "/features_norm.csv"))

# Dataframe with metainformation of contigs
coverage <- p$misc$ND_profile_path %>% 
  readRDS() %>% 
  attributes() %>% 
  `[[`("contig_coverage") %>% 
  left_join(
    get_GC_percent(p$misc$assembly_fasta), 
    by = c("chr" = "contig")) %>% 
  filter(chr %in% (features$contig %>% unique))

# Fasta sequences of contigs
fasta <- read.fasta(p$misc$assembly_fasta)

##-----------------------------------------------------------------------------
# GC vs. coverage plot
##-----------------------------------------------------------------------------
ggsave(
  filename = paste0(p$misc$output_folder, "/gc_vs_cov.pdf"),  
  width = 8, 
  height = 6,
  coverage %>% 
    ggplot(aes(x = avg_cov.dataset_B, y = gc, size = contig_length.x)) +
      geom_point(alpha = 0.2) +
      scale_x_log10(limits = c(4, 700)) +
      scale_size(
        range = c(1, 15), 
        limits = c(25000, 900000), 
        guide = guide_legend(title = "Contig Size")) +
      labs(x = "Mean Contig Coverage", y = "GC% of Contig") +
      theme_bw()
)

##-----------------------------------------------------------------------------
# t-SNE
##-----------------------------------------------------------------------------
source("./code/scripts/t-SNE.R")

##-----------------------------------------------------------------------------
# PCA
##-----------------------------------------------------------------------------
source("./code/scripts/PCA.R")

##-----------------------------------------------------------------------------
# Autoencoder
##-----------------------------------------------------------------------------

# system("
#   conda activate py38
#   python /shared-nfs/SH/code/scripts/AE.py
# ")

features_AE <- fread("/shared-nfs/SH/analysis_output/MGM1/test/binnning_AE_representation.tsv", sep = "\t")
features_AE <- features_AE[
  , cluster := hdbscan(.SD, minPts = 3)$cluster
][
  , cluster := ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
]
features_AE <- cbind(features_AE, tsne = features_AE[, Rtsne(.SD, theta = 0.3)$Y, .SDcols = 1:42]) %>% 
  cbind(chr = features$contig) %>% 
  merge(coverage, by = "chr")

# ggsave(
#   filename = paste0(p$misc$output_folder, "/cluster_AE.pdf"),  
#   width = 8, height = 6,
  ggplot(
    features_AE,
    aes(x = tsne.V1, y = tsne.V2, size = contig_length.x)) + 
    # geom_density_2d(
    #   contour = TRUE, 
    #   col = "gray60") +
    geom_point(
      aes(fill = cluster), 
      shape = 21, 
      alpha = 0.6) +
    theme_bw() +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2",
      title = "Autoencoder with 42 latent dimensions"
    ) +
    scale_size(
      range = c(0.1, 10)
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    guides(
      size = guide_legend(title = "Size"),
      fill = guide_legend(title = "Cluster")
    )
# )  
  #   width = 8, height = 6,
  ggplot(
    features_AE,
    aes(y = gc, x = avg_cov.dataset_B, size = contig_length.x)) + 
    # geom_density_2d(
    #   contour = TRUE, 
    #   col = "gray60") +
    geom_point(
      aes(fill = cluster), 
      shape = 21, 
      alpha = 0.6) +
    theme_bw() +
    labs(
      x = "Coverage",
      y = "GC%",
      title = "GC% vs. Coverage with methylation clusters"
    ) +
    scale_size(
      range = c(0.1, 10)
    ) +
    scale_x_log10(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    guides(
      size = guide_legend(title = "Size"),
      fill = guide_legend(title = "Cluster")
    )
  # )  
  # Saving cluster fasta
write_cluster_fasta(
  df = features_AE, 
  metagenome_fasta = fasta,
  path_out = paste0(p$misc$output_folder, "/bin_fasta_AE_PhD/")
)
  
# Run CheckM
system(
  paste0(
    "bash /shared-nfs/SH/code/scripts/cluster_quality_check.sh",
    " -i ", p$misc$output_folder, "/bin_fasta_AE_PhD",
    " -o ", p$misc$output_folder, "/bin_fasta_AE_PhD",
    " -x ", "fasta",
    " -t 40"
  )
) 
  
##-----------------------------------------------------------------------------
# AAE
##-----------------------------------------------------------------------------
read_delim("./analysis_output/MGM1/aae.tsv", col_names = FALSE) %>% 
  cbind(chr = features$contig) %>% 
  left_join(tsne_features) %>% 
  ggplot(aes( x = X1, y = X2, fill = cluster, size = contig_length.x)) +

  geom_point(
    aes(fill = cluster), 
    shape = 21, 
    alpha = 0.6) +
  theme_bw() +
  labs(
    x = "Latent Variable 1",
    y = "Latent Variable 2",
    title = "Advesarial Autoencoder"
  ) +
  scale_size(
    range = c(0.1, 15), 
    limits = c(1000, 1000000)
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


##-----------------------------------------------------------------------------
# Binning comparison
##-----------------------------------------------------------------------------

metabat_bins <- read.delim("nanodisco/MGM1/metabat/metabat", header = FALSE) %>% 
  mutate(V1 = as.factor(V1))
cluster_plot(
  ae_features %>% 
    left_join(metabat_bins, by = c("chr"="V1")) %>% 
    ggplot(aes(x = V1, y = V2.x, fill = V2.y, size = contig_length.x)),
  savepath = paste0(p$misc$output_folder, "/cluster_AE_metabat-color.pdf"),
  xlab = "Latent Variable 1", 
  ylab = "Latent Variable 2"
)

# Autoencoder plot with t-SNE cluster colors
cluster_plot(
  rbind(
    ae_features %>% 
      ungroup() %>% 
      mutate(
        cluster2 = tsne_features$cluster,
        cluster_type="t-SNE"),
    ae_features %>% 
      ungroup() %>% 
      mutate(
        cluster2 = cluster,
        cluster_type="Autoencoder")
  ) %>% 
    ggplot(aes(x = V1, y = V2, size = contig_length.x, fill = cluster2))+
    facet_wrap(~cluster_type),
  savepath = paste0(p$misc$output_folder, "/cluster_AE_tSNE-color.pdf"),
  dim = c(14, 6)
)


##-----------------------------------------------------------------------------
# Run checkM of cluster before proceeding
##-----------------------------------------------------------------------------
if (p$misc$run_checkm == "TRUE") {
  # Run CheckM
  system(
    paste0(
      "bash ./code/scripts/cluster_quality_check.sh",
      " -i ", p$misc$output_folder, "/t_SNE_cluster_fasta",
      " -o ", p$misc$output_folder, "/t_SNE_cluster_fasta",
      " -x ", "fasta",
      " -t 20"
    )
  ) 
  
  
  # CheckM results of clusters 
  cluster_plot_checkm  <- read_checkm(paste0(p$misc$output_folder, "/t_SNE_cluster_fasta/storage/bin_stats_ext.tsv")) %>% 
    select(bin, Completeness, Contamination, Genome_size) %>% 
    mutate(
      Completeness = round(as.numeric(Completeness, digits = 4)),
      Contamination = round(as.numeric(Contamination, digits = 4)),
      comp_cont = paste0(str_pad(Completeness, 3, side = "right"), "; ", format(Contamination, 3, side = "right")),
      cluster = str_extract(bin, "(?<=\\_).*")
    ) %>% 
    select(-bin) %>% 
    mutate(cluster_2 = ifelse(as.numeric(Genome_size) > 5e1, paste0(
      str_pad(cluster, 2, side = "right"), 
      " - " ,
      comp_cont), NA)) %>% 
    right_join(tsne_features) %>% 
    ggplot(aes(x=V1, y=V2, fill=cluster_2, size = contig_length.x, fill = cluster_2)) +
      geom_point( 
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
    ggsave(cluster_plot_checkm, filename = paste0(p$misc$output_folder, "/cluster_tsne_checkm.pdf"))
    
    gtdbtk <- read_delim("analysis_output/MGM1/t_SNE_cluster_fasta/GTDBTk/gtdbtk.bac120.summary.tsv") %>% 
      group_by(user_genome) %>% 
      mutate(
        classification %>% 
          str_split(pattern = ";") %>% 
          unlist() %>% 
          str_remove(".\\_\\_") %>% 
          t() %>% 
          data.frame(stringsAsFactors = FALSE) %>% 
          setNames(c("d", "p", "c", "o", "f", "g", "s"))
      )
    #"/shared-nfs/SH/analysis_output/MGM1/bin_fasta_AE_PhD/storage/bin_stats_ext.tsv"
    a <- read_checkm(paste0(p$misc$output_folder, "/t_SNE_cluster_fasta/storage/bin_stats_ext.tsv")) %>% 
      select(bin, Completeness, Contamination, Genome_size, contigs) %>% 
      mutate(
        Completeness = round(as.numeric(Completeness, digits = 4)),
        Contamination = round(as.numeric(Contamination, digits = 4)),
        Genome_size = paste0(round(as.numeric(Genome_size)/10^6, 1), " Mb"),
        bin = str_extract(bin, "(?<=\\_).*")
      ) %>% 
      arrange(as.numeric(bin)) %>% 
      gt::gt() %>% 
      gt::data_color(
        columns = Completeness,
        colors = scales::col_numeric(
          palette = c("red2", "green2"),
          domain = c(0,100)),
        alpha = 0.3,
        autocolor_text = FALSE
      ) %>% 
      gt::data_color(
        columns = Contamination,
        colors = scales::col_numeric(
          palette = c("green2", "red2"),
          domain = c(0,20), na.color = "red2"),
        alpha = 0.3, 
        autocolor_text = FALSE
      ) %>% 
      gt::data_color(
        columns = bin,
        alpha = 0.8,
        colors =c(
          '#e6194B', #1
          'gray50',  #10
          '#aaffc3', #11
          '#3cb44b', #2
          '#ffe119', #3
          '#4363d8', #4
          '#f58231', #5
          '#fabed4', #6
          '#469990', #7
          '#dcbeff', #8
          '#fffac8'  #9
          )
      ) %>% 
      gt::cols_align(
        columns = everything(),
        align = "center"
      ) %>% 
      gt::cols_align(
        columns = bin,
        align = "left"
      ) %>% 
      gt::cols_label(
        bin = "Cluster", 
        Genome_size = "Genome Size", 
        contigs = "Nr. Contigs"
        ) %>% 
      gt::tab_style(
        style =list(
          gt::cell_borders(
            sides = c("right", "left", "top", "bottom"),
            color = "#606060",
            weight = gt::px(2)
          )
        ),
        locations = list(gt::cells_body(
          columns = everything()
        ))
      )
    pdf(paste0(p$misc$output_folder, "/checkm_table.pdf"), a)
    a
    ggsave(
      filename = paste0(p$misc$output_folder, "/checkm_table.pdf"),  
      width = 8, height = 6,
      a)
    }

##-----------------------------------------------------------------------------
# Investigating the effect of features 
##-----------------------------------------------------------------------------

# Extracting features with higher within cluster mean and standard deviation
high_info_features <- features %>% 
  left_join(
    select(tsne_features, chr, cluster), 
    by = c("contig"="chr")) %>% 
  pivot_longer(
    cols = c(-contig, -cluster), 
    names_to = "feature", 
    values_to = "dist_score") %>% 
  group_by(feature, cluster) %>% 
  summarise(
    variance = var(dist_score),
    mean = mean(abs(dist_score)),
    sd = sd(dist_score),
    z_score = sd/mean) %>% 
  left_join(
    features %>% 
      left_join(
        select(tsne_features, chr, cluster), 
        by = c("contig"="chr")) %>% 
      pivot_longer(
        cols = c(-contig, -cluster), 
        names_to = "feature", 
        values_to = "dist_score") %>% 
      group_by(feature) %>% 
      summarise(
        variance_all = var(dist_score),
        mean_all = abs(mean(dist_score)),
        sd_all = sd(dist_score),
        z_score_all = sd_all/mean_all)) %>% 
  mutate(
    sd_dif = sd - sd_all,
    variance_dif = variance - variance_all,
    mean_dif = mean - mean_all,
    z_score_dif = z_score_all
  ) %>% 
  #ilter(mean_dif > 1) %>% 
  mutate(motif = str_extract(feature, ".*(?=\\_)"))

# Complete range of difference
ggsave(
  filename = paste0(p$misc$output_folder, "mean-dif_vs_var-dif_full_range.pdf"),  
  width = 12, height = 8,
  ggplot(
      high_info_features, 
      aes(x = variance_dif, y = mean_dif)) +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = cluster), size = 1, alpha = 0.5) +
    #geom_text(aes(label = feature), size = 1.5, position = "dodge", angle = 15) +
    facet_wrap(~cluster, scale = "free") + 
    labs(
      x = "Cluster Variance - Total variance",
      y = "Cluster Mean - Total Mean") +
    theme_bw()
)

# Selected range of differences
ggsave(
  filename = paste0(p$misc$output_folder, "mean-dif_vs_var-dif_zoomed.pdf"),  
  width = 15, height = 8,
  ggplot(high_info_features %>% 
           group_by(cluster) %>% 
           arrange(-mean_dif) %>% 
           filter( mean_dif > 1 & sd_dif < 0), 
         aes(x = variance_dif, y = mean_dif)) +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 0) +
    #geom_point(aes(col = cluster), size = 1, alpha = 0.5) +
    geom_text(aes(label = feature, col = motif), size = 2) +
    facet_wrap(~cluster) + 
    labs(
      x = "Cluster Variance - Total variance",
      y = "Cluster Mean - Total Mean") +
    theme_bw()
)

# Checking if removal of these feature directly remove cluster in tSNE
dup_row <- features %>% 
  `[`(!(names(features) %in% (high_info_features %>% filter(mean_dif > 1 & cluster == "5") %>% pull(feature) %>% unique()))) %>% 
  select(-contig) %>%
  duplicated() %>% 
  which()
features %>% 
  `[`(!(names(features) %in% "AGATG_1")) %>% 
  select(-contig) %>%
  distinct() %>% 
  tsne() %>% 
  mutate(
    cluster = hdbscan(cbind(V1, V2), minPts = 18)$cluster,
    cluster = ifelse(cluster == 0, NA, cluster) %>% factor(levels = sort(unique(.)))
  ) %>%
  cbind(chr = features$contig) %>% 
  left_join(coverage) %>% 
  #left_join(tsne_checkm, by = "chr") %>% 
  ggplot(aes(x = V1, y = V2, size = contig_length.x)) + 
    stat_density_2d(col = "gray60") +
    geom_point(aes(fill = cluster), shape = 21, alpha = 0.9) +
    scale_size(range = c(1, 15), limits = c(25000, 900000)) +
    labs(x = "t-SNE dim 1", y = "t-SNE dim 2", title = "") +
    theme_bw()





# Evaluating feature shared by clusters
high_info_features %>% 
  group_by(cluster) %>% 
  arrange(-mean_dif) %>% 
  filter( mean_dif > 1 & sd_dif < 0) %>%
  pull(feature) %>% 
  unique() %>% 
  `[`(features, .) %>% 
  cbind(chr = coverage$chr) %>% 
  pivot_longer(-chr, names_to = "feature", values_to = "dist_score") %>% 
  left_join(tsne_features, by = "chr") %>% 
  ggplot(aes(x = feature, y = dist_score)) +
  #geom_jitter(aes(col = cluster), alpha = 0.5) +
  geom_boxplot(
    width = 0.7,
    alpha = 1,
    outlier.size = 0.6,
    outlier.color = NA
  ) + 
  facet_wrap(~cluster, scales = "free") +
  labs(y = "Distance Score", x = "Motif_Distance", title = "Motif mean dif > 1 and sd dif < 0") +
  coord_flip() +
  theme_bw()


# Looking at distrubution of a feature for each cluster
dir.create(
  paste0(p$misc$output_folder, "/cluster_distrubution_of_motifs"),
  showWarnings = FALSE
  )
high_info_features %>% 
  filter(mean_dif > 1 & sd_dif < 0) %>% 
  pull(feature) %>% 
  unique() %>% 
  lapply(function(feature) {
    plot <- features %>% 
      `[`(feature) %>% 
      cbind(chr = coverage$chr) %>% 
      pivot_longer(-chr, names_to = "feature", values_to = "dist_score") %>% 
      left_join(tsne_features, by = "chr") %>% 
      filter(contig_length.x > 10000 & bin_size > 1e6) %>% 
      ggplot(
          aes(
            x = cluster, 
            y = dist_score, 
            fill = paste0(
                    "Cluster ", cluster, 
                    "\nSize: ", format(bin_size, scientific = TRUE, big.mark = ","),
                    "\n", n_contigs)
            )
          ) +
        ggdist::stat_halfeye(
          .width = 0,
          point_color = NA,
          justification = -0.25
          ) +
        geom_jitter(
          shape = 21, 
          width = 0.1, 
          alpha = 0.7, 
          size = 0.9
          ) +
        geom_boxplot(
          width = 0.5, 
          outlier.colour = NA
          ) +
        labs(
          x = "Cluster", 
          y = "Distance score", 
          title = paste0("Distrubution of motif; ", feature)
          ) +
        coord_flip() +
        theme_bw()
    ggsave(
      plot,
      filename = paste0(p$misc$output_folder, "/cluster_distrubution_of_motifs/", feature, ".pdf"),
      device = "pdf",
      width = 15, height = 8
    )
})

