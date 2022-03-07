#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/software/Rlib" , .libPaths()))

setwd(here::here())
pacman::p_load(
  "data.table",
  "dplyr",
  "ggplot2", 
  "rhdf5",
  "seqinr",
  "Rsamtools",
  "stringr",
  "patchwork",
  "logger",
  "glue",
  "pbapply",
  "scales"
)

###############################################################################
# Utility functions
###############################################################################
# goup_nest version for data.table
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

# unnest version for data.table
unnest_dt <- function(dt, col, id){
  stopifnot(is.data.table(dt))
  by <- substitute(id)
  col <- substitute(unlist(col, recursive = FALSE))
  dt[, eval(col), by = eval(by)]
}

set_sum_1 <- function(x){
  x/sum(x)
}

# Calculate rolling mean with distance weigth decay
rolling_mean <- function(x, n = 2, weight_dropoff = 0.75){
  stopifnot(is.numeric(x))
  max_i <- length(x)
  weigths <- lapply(1:n, function(x) weight_dropoff^x) %>% unlist()
  weigths <- c(weigths[length(weigths):1], 1, weigths)
  out <- lapply(
      seq_along(x),
      function(i){
        values <- as.numeric(x[max(1,i-n):min(max_i, i+n)])
        weigths[max(1, n - i + 2):min(n*2 + 1, n + 1 + max_i - i)] %>% 
          set_sum_1() %>%  
        sum(values * ., na.rm = TRUE)
      }
    ) %>%  unlist()
  return(out)  
}

# Calculate rolling mean in log space with distance weigth decay
rolling_mean_log <- function(x, contig_index, n = 2, weight_dropoff = 0.75, base = 10, data_splits = 25){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(contig_index))
  warningCondition(length(x) < 1, message = "Input vector of size 0?")
  warningCondition(length(contig_index) < 1, message = "Index vector of size 0?")
  warningCondition(is.na(contig_index) %>%  any(), message = "There are NA values in the contig index")

  dir.create(paste0(arg$out, "/temp"), showWarnings = FALSE)
  if(length(na.omit(x)) != 0) {
    log_info(glue("length of x: {length(x)}"))
    log_info(glue("length of x, no NA: {length(na.omit(x))}"))
    log_info(glue("Percent NA: {1 - (length(na.omit(x)) / length(x))}"))
    log_info(glue("Mean P-values: {mean(x, na.rm = TRUE)}"))
    pboptions(nout = data_splits)
    split_ind <- c(rep(1:data_splits, each = floor(length(x)/data_splits)), 
      rep(data_splits, length(x) %% data_splits))
    x_batches <- split(x, split_ind)
    pos_batches <- split(contig_index, split_ind)
    for (i in 1:data_splits) {
      chunk_x <- x_batches[[i]]
      save(chunk_x, file = paste0(arg$out, "/temp/x_bacth_", i, ".Rdata"))
      chunk_pos <- pos_batches[[i]]
      save(chunk_pos, file = paste0(arg$out, "/temp/pos_bacth_", i, ".Rdata"))
    }
    max_ind <- max(contig_index)
    weigths <- lapply(1:n, function(x) weight_dropoff^x) %>% unlist()
    weigths <- c(weigths[length(weigths):1], 1, weigths)
    out <- pblapply(
      1:data_splits,
      function(i) {
        load(file = paste0(arg$out, "/temp/x_bacth_", i, ".Rdata"))
        load(file = paste0(arg$out, "/temp/pos_bacth_", i, ".Rdata"))
        lapply(
          chunk_pos,
          function(pos){
            ind <- max(1, pos - n):min(max_ind, pos + n)
            weigths[chunk_pos[(chunk_pos %in% ind)] - pos + n + 1] %>% 
              set_sum_1() %>%  
              sum(log(chunk_x[chunk_pos %in% ind], base = base) * ., na.rm = TRUE) %>% 
              `^`(10, .)
          }
        )
      },
      cl = arg$threads
    ) %>%  unlist()
    log_info(glue("Mean weighted P-values: {mean(out)}"))
    return(out)
  } else {
    log_info("All values in vector are NA")
    return(rep(NA, length(x)))
  }
}

###############################################################################
# Input arguments
###############################################################################

arg <- list()
arg$reference <- "/shared-nfs/SH/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta"
arg$signal_mapping_pcr <- "/shared-nfs/SH/samples/zymo/megalodon/pcr_test/signal_mappings.hdf5"
arg$signal_mapping_nat <- "/shared-nfs/SH/samples/zymo/megalodon/nat/signal_mappings.hdf5"
arg$read_mapping_pcr <- "/shared-nfs/SH/samples/zymo/megalodon/pcr_test/mappings.sort.view.txt"
arg$read_mapping_nat <- "/shared-nfs/SH/samples/zymo/megalodon/nat/mappings.sort.view.txt"
arg$out <- "/shared-nfs/SH/samples/zymo/current_difference/two_contigs"
arg$contig_plot <- c("bs_contig1", "lf_contig1")
arg$overwrite <- FALSE
arg$save_intermediate_files <- FALSE
arg$min_u_val <- 1e-12
arg$u_val_weight_window <- 3 # val*2 + 1
arg$u_val_weight_dropoff <- 0.9
arg$min_mappings <- 3
arg$threads <- 20

dir.create(arg$out, showWarnings = FALSE)
jsonlite::write_json(arg, path = glue("{arg$out}/args.json"), pretty = TRUE)



###############################################################################
##                                                                           ##
##                                                                           ##
##                     Functions                                             ##
##                                                                           ##
##                                                                           ##
###############################################################################

load_mapping_hdf5 <- function(hdf5_file, 
                              batch,
                              mappings,
                              type,
                              reads_ids = NA){
  # Read dacs (current measurements)
  log_info("Loading data from {batch}")
  dacs <- h5read(hdf5_file, name = paste0("/Batches/", batch, "/Dacs")) %>% 
    data.table() %>% 
    setnames("dacs")
  
  # Read Ref_to_signal (Number of dacs mapping to each ref position)
  Ref_to_signal <- h5read(hdf5_file, name = paste0("/Batches/", batch, "/Ref_to_signal")) %>% 
    data.table() %>% 
    setnames("ref_to_signal")
  
  # Read Reference (Interger representation of read bases)
  Reference = h5read(hdf5_file, name = paste0("/Batches/", batch, "/Reference")) %>% 
    data.table() %>% 
    setnames("reference")
  
  # Read read information (mainly ID, dacs length, Ref/sig length)
  reads <- h5read(hdf5_file, name = paste0("/Batches/", batch, "/Dacs_lengths")) %>% 
    data.table() %>% 
    setnames("Dacs_lengths") %>% 
    cbind("Ref_to_signal_lengths" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/Ref_to_signal_lengths"))) %>% 
    cbind("Reference_lengths" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/Reference_lengths"))) %>% 
    #cbind("digitisation" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/digitisation"))) %>% 
    cbind("offset" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/offset"))) %>% 
    cbind("range" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/range"))) %>% 
    #cbind("scale_frompA" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/scale_frompA"))) %>% 
    #cbind("shift_frompA" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/shift_frompA"))) %>% 
    apply(MARGIN = 2, FUN = as.numeric) %>% 
    data.table(stringsAsFactors = FALSE) %>% 
    cbind("read_id" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/read_id")))

  
  log_info("Adding read information to reference mapping")
  # Add read id to ref to signal for later grouping
  Ref_to_signal <- mapply(
      function(len, start, id){
        list(
          signal_frame = Ref_to_signal$ref_to_signal[start:(start+len-1)],
          read_id = id
        )
      },
      reads$Ref_to_signal_lengths,
      (cumsum(reads$Ref_to_signal_lengths) - reads$Ref_to_signal_lengths + 1),
      reads$read_id, 
      SIMPLIFY =  FALSE
    ) %>% 
      rbindlist()
  
  # Add dacs window and reference integer
  Ref_to_signal[
      , n_dac := signal_frame-data.table::shift(signal_frame, fill = 0), by = read_id
    ]
  Ref_to_signal <- na.omit(Ref_to_signal, cols=seq_along(Ref_to_signal))[
      , start := cumsum(n_dac)
    ][
      !(n_dac == signal_frame & n_dac != 0),
    ][
      , ref := Reference$reference
    ][
      , pos := 1:.N, by = read_id
    ]
  log_info("Adding reference mapping to dacs")
  
  if (!is.na(reads_ids)) {
    Ref_to_signal <- Ref_to_signal %>% 
      subset(read_id %in% reads_ids)
  }
  dacs <- dacs$dacs
  Ref_to_signal <- Ref_to_signal[
    , dac := mapply(
      function(len, end){
        dac <- dacs[(end-len+1):end]
      },
      n_dac,
      start,
      SIMPLIFY = FALSE
    ) 
  ][
    , batch := batch
  ][
    , unlist(dac, recursive = FALSE), 
    by = eval(names(Ref_to_signal)[!(names(Ref_to_signal) %in% c("dac"))])
  ] %>% 
    setnames("V1", "dac")
  
  add_mapping_to_dacs(Ref_to_signal, mappings = mappings, type = type)
  
  log_success("Finished processing {batch}")
  return(Ref_to_signal)
}


###############################################################################
# Process signal mappings
###############################################################################

# Append reference mappings
add_mapping_to_dacs <- function(dacs, mappings, type){
  log_info("Joining mappings")
  dacs[
      mappings, on = .(read_id), `:=`(contig = i.contig, start = i.start, direction = i.direction)
    ][
      , type := type
    ][
      , pos := pos - 1
    ][
      , contig_index := start + pos
    ][
      , contig_id := paste0(contig, "_", contig_index)
    ]
  log_info("Normalising dacs")
  dacs[
      , dacs_norm := (dac - mean(dac))/sd(dac), by = read_id
    ]
}

# Calculate mean difference and u-test value
calculate_current_diff <- function(
  dt1, 
  dt2, 
  min_u_val, 
  min_cov = 3, 
  data_splits = 1,
  parallel = TRUE){
  # log_info("Setting up cluster")
  # cl <- makeForkCluster(nnodes = arg$threads)
  pboptions(nout = data_splits)
  log_info("Mean difference")
  dt1[
      dt2, on = "contig_id", nat := i.nat
    ][
      , mean_nat := pblapply(
          nat,          
          function(x){
            if (is.null(x)) {
              NA
            } else {
              mean(x$dacs_norm, na_rm = TRUE)}
          },
          cl = ifelse(parallel, arg$threads, NULL)) %>% 
        unlist()
    ][
      , mean_pcr := pblapply(
          pcr, 
          function(x){
            if (is.null(x)) {
              NA
            } else {
              mean(x$dacs_norm, na_rm = TRUE)}
            },
          cl = ifelse(parallel, arg$threads, NULL)) %>% 
        unlist()
    ][
      , mean_dif := mean_nat - mean_pcr
    ]
  log_info("Dac number and coverage")
  parallel = FALSE
  dt1[
      , n_nat := pblapply(
          nat, 
          nrow,
          cl = if(parallel) arg$threads else NULL) %>% 
        pblapply(
          X = .,
          function(x){ifelse(is.null(x), NA, x)},
          cl = if(parallel) arg$threads else NULL) %>% 
        unlist()
    ][
      , n_pcr := pblapply(
          pcr,
          nrow,
          cl = if(parallel) arg$threads else NULL) %>% 
        pblapply(
          X = .,
          function(x){ifelse(is.null(x), NA, x)},
          cl = if(parallel) arg$threads else NULL) %>% 
        unlist()
    ][
      , n_nat_map := pblapply(
          nat, 
          function(x){x$read_id %>% unique %>% length},
          cl = if(parallel) arg$threads else NULL) %>% 
        unlist()
    ][
      , n_pcr_map := pblapply(
          pcr, 
          function(x){x$read_id %>% unique %>% length},
          cl = if(parallel) arg$threads else NULL) %>% 
        unlist()
    ]
  parallel = TRUE
  log_info("Wilcox test")
  dt1[
      , u_val := mcmapply(
        function(nat, pcr, n_nat_map, n_pcr_map){
            if(n_nat_map >= min_cov & n_pcr_map >= min_cov) {
              wilcox.test(nat$dacs_norm, pcr$dacs_norm)$p.value
            } else {
              NA
            }
          },
        nat,
        pcr,
        n_nat_map,
        n_pcr_map,
        mc.cores = arg$threads
      ) %>% unlist()
    ]
  log_info("Weighting p-values")
  dt1[
      , u_val_weighted_log := rolling_mean_log(u_val, contig_index = contig_index, n = arg$u_val_weight_window, weight_dropoff = arg$u_val_weight_dropoff), by = contig
    ][
      , u_val_bonf := p.adjust(u_val, n = length(u_val), method = "bonferroni")
    ][
      , u_val_bonf_weighted_log := rolling_mean_log(u_val_bonf, contig_index = contig_index, n = arg$u_val_weight_window, weight_dropoff = arg$u_val_weight_dropoff), by = contig
    ][
      , u_val_BH := p.adjust(u_val, n = length(u_val), method = "BH")
    ][
      , u_val_BH_weighted_log := rolling_mean_log(u_val_BH, contig_index = contig_index, n = arg$u_val_weight_window, weight_dropoff = arg$u_val_weight_dropoff), by = contig
    ]
}

###############################################################################
# Visualisation of events in one contig
###############################################################################
# Plots 
pboptions(type = "timer")
plot_events <- function(dt){
  guides_format <- guides(
      fill = guide_legend(title = "Event"),
      size = guide_legend(title = "Event"),
      col = guide_legend(title = "Event"),
      alpha = "none"
    )
  plot_col <- c("#99712d", "#2d5599")
  geom_point_default <- geom_jitter(alpha = 0.1, size = 0.1)
  
  
  p1.1 <- dt %>% 
    ggplot(aes(x = contig_index, y = mean_dif)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = mean_dif), y = 0, size = 0.2) +
    geom_point(shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (Mean difference)"
    )  + 
    scale_size_manual(values = c(0.4, 2)) +
    scale_fill_manual(values = plot_col) +
    guides_format +
    ggtitle("mean_diff")
  
  p1.2 <- dt %>% 
    ggplot(aes(x = contig_index, y = u_val_weighted_log)) +
    geom_hline(yintercept = 1) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = u_val_weighted_log), y = 0, size = 0.2) +
    geom_point(shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (U-value)"
    ) +
    scale_size_manual(values = c(0.4, 2)) +
    scale_y_log10() +
    scale_fill_manual(values = plot_col) +
    guides_format +
    ggtitle("p_val")
    
  
  p1.3 <- dt %>%  
    ggplot(aes(x = n_pcr, y = n_nat)) +
    #geom_point_default +
    labs(
      x = "PCR dacs",
      y = "NAT dacs"
    ) + 
    geom_point_default +
    stat_density_2d_filled(
      contour = TRUE,
      aes(alpha = after_stat(level)),
      bins = 10, col = "gray10", size = 0.1
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_alpha_discrete(range = c(0,1)) +
    scale_fill_manual(values = plot_col) +
    scale_color_manual(values = plot_col) +
    guides_format +
    ggtitle("dacs")
  
  p1.4 <- dt %>% 
    ggplot(aes(x = n_pcr_map, y = n_nat_map, col = event)) +
    labs(
      x = "PCR coverage",
      y = "NAT coverage"
    ) + 
    geom_point_default +
    # stat_density_2d_filled(
    #   contour = TRUE,
    #   aes(alpha = after_stat(level)),
    #   bins = 2, col = "gray10", size = 0.1
    # ) +
    scale_alpha_discrete(range = c(0,1)) +
    scale_color_manual(values = plot_col) +
    scale_fill_manual(values = plot_col) +
    guides_format +
    ggtitle("coverage")
  
  p1.5 <- dt %>% 
    mutate(
      mappings = n_pcr_map / n_nat_map
    ) %>% 
    ggplot(aes(y = mappings, x = u_val_weighted_log)) +
    labs(
      x = "U-value",
      y = "nr. PCR / nr. NAT"
    ) + 
    geom_point_default +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format +
    ggtitle("cov_proportion")
  
  p1.6 <- dt %>% 
    ggplot(aes(x = mean_dif, y = u_val_weighted_log)) +
    labs(
      x = "Mean difference",
      y = "U-value"
    ) + 
    geom_point_default +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format +
    ggtitle("fake_vulcano_plot")

  return(list(p1.1, p1.2, p1.3, p1.4, p1.5, p1.6))
}

###############################################################################
##                                                                           ##
##                                                                           ##
##                     Processing                                            ##
##                                                                           ##
##                                                                           ##
###############################################################################

# Importing reads to reference mappings
mapping_pcr <- fread(arg$read_mapping_pcr) %>% 
  select(V1, V3, V4, V2) %>% 
  setnames(c("read_id", "contig", "start", "direction")) %>% 
  mutate(
    direction = if_else(direction == 16, "fwd", "rev")
  )

mapping_nat <- fread(arg$read_mapping_nat) %>% 
  select(V1, V3, V4, V2) %>% 
  setnames(c("read_id", "contig", "start", "direction")) %>% 
  mutate(
    direction = if_else(direction == 16, "fwd", "rev")
  )

###############################################################################
# PCR reads
###############################################################################

if (TRUE) {
  # Loading signal mappings
  log_info("Processing PCR mappings")
  hdf5_pcr <- H5Fopen(arg$signal_mapping_pcr)
  signal_mappings_pcr_unnested <- h5ls(hdf5_pcr) %>% filter(group == "/Batches") %>% pull(name) %>% 
    pblapply(
      mc.preschedule = FALSE,
      cl = arg$threads,
      function(batch){
        load_mapping_hdf5(
          hdf5_pcr, 
          batch = batch,
          mappings = mapping_pcr, 
          type = "pcr",
          reads_ids = mapping_pcr[contig %in% arg$contig_plot,][, read_id]
        )
      }
    ) %>% 
    rbindlist()
  h5closeAll()
  rm(hdf5_pcr)

  if (arg$save_intermediate_files) {
    # Save signal mappings
    log_info("Saving processed PCR dacs")
    fwrite(signal_mappings_pcr_unnested, glue("{arg$out}/singal_mappings_pcr.tsv"))
  }
  
  signal_mappings_pcr <- signal_mappings_pcr_unnested %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "pcr")
  log_success("Succesfully processed PCR dacs")

} 

###############################################################################
# NAT reads
###############################################################################

if (TRUE) {
  log_info("Processing NAT signal mappings")
  hdf5_nat <- H5Fopen(arg$signal_mapping_nat)
  signal_mappings_nat_unnested <- h5ls(hdf5_nat) %>% filter(group == "/Batches") %>% pull(name) %>% 
    pblapply(
      mc.preschedule = FALSE,
      cl = arg$threads,
      mc.silent = FALSE,
      function(batch){
        load_mapping_hdf5(
          hdf5_nat, 
          batch = batch,
          mappings = mapping_nat, 
          type = "nat",
          reads_ids = mapping_nat[contig %in% arg$contig_plot,][, read_id]
        )
      }
    ) %>% 
    rbindlist()
  h5closeAll()
  rm(hdf5_nat)

  if (arg$save_intermediate_files) {
    # Save signal mappings
    log_info("Saving processed NAT dacs")
    fwrite(signal_mappings_nat_unnested, glue("{arg$out}/singal_mappings_nat.tsv"))
  }
  signal_mappings_nat <- signal_mappings_nat_unnested %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "nat")
  log_success("Succesfully processed NAT dacs")
  
}

###############################################################################
# Calculate current difference
###############################################################################

if (!file.access(glue("{arg$out}/current_difference.tsv"), mode = 4) == 0 | arg$overwrite) {
  log_info("Calculating current difference and statistics")
  current_difference <- copy(signal_mappings_pcr)
  calculate_current_diff(
      current_difference, 
      signal_mappings_nat, 
      arg$min_u_val
    )
  current_difference[, `:=`(nat = NULL, pcr = NULL)]
   
  # Save processed signal mappings
  log_info("Saving current difference")
  fwrite(
    current_difference,
    glue("{arg$out}/current_difference.tsv")
  )
  log_success("Performed currrent difference")
} else {
  current_difference <- fread(glue("{arg$out}/current_difference.tsv"))
}

###############################################################################
# Plot
###############################################################################

# p_all <- plot_events(current_difference)
# 
# p_all %>% 
# pblapply(
#   function(p){
#       ggsave(
#         width = unit(8, "cm"),
#         height = unit(10, "cm"),
#         p, 
#         filename = glue("{arg$out}/current_difference_{p$labels$title}.png"), 
#         device = "png"
#       )
#     }
#   )


###############################################################################
# Evaluation of PCR batch 1 vs batch 0
###############################################################################

signal_mappings_pcr_0 <- signal_mappings_pcr_unnested[batch == "Batch_0", ] %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "pcr")

signal_mappings_pcr_1 <- signal_mappings_pcr_unnested[batch == "Batch_1", ] %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "nat")

current_difference_pcr <- copy(signal_mappings_pcr_0)
calculate_current_diff(current_difference_pcr, signal_mappings_pcr_1, arg$min_u_val)


# plot_events(signal_mappings_pcr_0[
#   n_nat_map > 1 & n_pcr_map > 1
# ])

# ggsave(
#   width = unit(8, "cm"),
#   height = unit(10, "cm"),
#   p_pcr_vs_pcr, 
#   filename = glue("{arg$out}/current_difference_pcr_vs_pcr_{arg$contig_plot}.png"), 
#   device = "png"
# )

###############################################################################
# Evaluation of differenct U-value cut offs
###############################################################################

# pcr vs. pcr and nat vs. pcr plot function
plot_PvP_NvP <- function(thresholds, col){
  thresholds %>% 
    lapply(
      function(x){
        n_NvP <- sum(abs(current_difference[, ..col]) < x, na.rm = TRUE)
        p_NvP <- n_NvP / nrow(na.omit(current_difference[, ..col]))
        n_PvP <- sum(abs(current_difference_pcr[, ..col]) < x, na.rm = TRUE)
        p_PvP <- n_PvP / nrow(na.omit(current_difference_pcr[, ..col]))
        list(
          x,
          p_NvP, 
          p_PvP,
          p_NvP - p_PvP,
          p_PvP/p_NvP
        )
      }) %>% 
    rbindlist() %>% 
    setnames(paste0("V", 1:5), c("U_val", "p_NvP",  "p_PvP", "p_diff", "p_ratio")) %>%  
    tidyr::pivot_longer(starts_with("p_"), values_to = "proportion", names_to = "type") %>% 
    ggplot() +
    aes(x = U_val, y = proportion, col = type) +
    geom_line(size = 1.2) +
    #geom_point(shape = 21, size = 1.5) +
    scale_x_log10() +
    scale_y_continuous(labels = percent) +
    labs(
      x = "P-value threshold",
      y = "Proportion of events"
    ) +
    theme_bw() +
    scale_color_manual(values = c("#127474", "#741212", "#437412", "#431274"))
}
ylimits <- c(1e-5, 1)
ybreaks <- 10 ^ seq(log(ylimits[2], base = 10), log(ylimits[1], base = 10), by = -1)
ylabel <- paste0(signif(ybreaks*100, 1), "%")
# p-value
p2.1 <-  10^-seq(0, 20, by = 0.2) %>% 
  plot_PvP_NvP(col = "u_val") +
  labs(title = "P-value") +
  scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

# p-value log weighted 
p2.3 <-  10^-seq(1, 20, by = 0.2) %>% 
  plot_PvP_NvP(col = "u_val_weighted_log") +
  labs(title = "P-value log weighted") +
  scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)
  
# Mean difference cutoff
p2.4 <- 10^-seq(-1, 3, by = 0.1) %>% 
  plot_PvP_NvP("mean_dif") +
  xlab("Mean difference threshold") +
  labs(title = "Mean difference")

p2.5 <- 10^-seq(0, 20, by = 0.1) %>% 
  plot_PvP_NvP("u_val_BH") +
  labs(title = "P-value B&H corrected") +
  scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

p2.6 <- 10^-seq(0, 20, by = 0.1) %>% 
  plot_PvP_NvP("u_val_BH_weighted_log") +
  labs(title = "P-value B&H corrected and log weighted") +
  scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)

p2.7 <- 10^-seq(0, 20, by = 0.1) %>% 
  plot_PvP_NvP("u_val_bonf") +
  labs(title = "P-value Bonferroni") +
  scale_y_log10(breaks = ybreaks, limits = ylimits, label = ylabel)


p2_all <- (p2.1 + p2.5 + p2.7) /
(p2.3 + p2.6 + p2.4) +
plot_layout(guides = "collect") +
  plot_annotation(
    title =  "Estimation of event cut-off"
  ) & 
  theme_bw()
ggsave(
  width = unit(16, "cm"),
  height = unit(7, "cm"),
  filename = glue("{arg$out}/event_cutoff_estimation.png"),
  p2_all,
  device = "png")





