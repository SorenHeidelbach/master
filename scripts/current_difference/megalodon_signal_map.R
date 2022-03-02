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
  "pbapply"
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

# Calculate rolling mean with distance weigth decay
rolling_mean <- function(x, n = 2, weigth_dropoff = 0.75){
  max_i <- length(x)
  weigths <- lapply(1:n, function(x) weigth_dropoff^x) %>% unlist()
  weigths <- c(weigths[length(weigths):1], 1, weigths)
  out <- c()
    for (i in seq_along(x)) {
      values <- x[max(1,i-n):min(max_i, i+n)]
      weigths_indexed <- weigths[max(1, n - i + 2):min(n*2 + 1, n + 1 + max_i - i)]
      weigths_indexed <- weigths_indexed/sum(weigths_indexed)
      out <- c(out, sum(values * weigths_indexed))
    }
  return(out)  
}

# Calculate rolling mean in log space with distance weigth decay
rolling_mean_log <- function(x, n = 2, weigth_dropoff = 0.75, base = 10){
  max_i <- length(x)
  weigths <- lapply(1:n, function(x) weigth_dropoff^x) %>% unlist()
  weigths <- c(weigths[length(weigths):1], 1, weigths)
  out <- c()
  for (i in seq_along(x)) {
    values <- log(x[max(1,i-n):min(max_i, i+n)], base = base)
    weigths_indexed <- weigths[max(1, n - i + 2):min(n*2 + 1, n + 1 + max_i - i)]
    weigths_indexed <- weigths_indexed/sum(weigths_indexed)
    out <- c(out, sum(values * weigths_indexed))
  }
  out <-base ^ out
  return(out)  
}

###############################################################################
# Input arguments
###############################################################################

arg <- list()
arg$reference <- "/shared-nfs/SH/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta"
arg$signal_mapping_pcr <- "/shared-nfs/SH/samples/zymo/megalodon/pcr_test/signal_mappings.hdf5"
arg$signal_mapping_nat <- "/shared-nfs/SH/samples/zymo/megalodon/nat_test2/signal_mappings.hdf5"
arg$out <- "/shared-nfs/SH/samples/zymo/current_difference"
arg$contig_plot <- "lf_contig3"
arg$min_u_val <- 1e-12
arg$u_val_weight_window <- 4 # val*2 + 1
arg$u_val_weight_dropoff <- 0.9
arg$overwrite <- TRUE
arg$min_mappings <- 3

###############################################################################
# Import signal mappings
###############################################################################

# Load in reference
ref <- read.fasta(arg$reference)

load_mapping_hdf5 <- function(hdf5_file, batch, reads_ids = NA){
  # Read dacs (current measurements)
  log_info("Loading data from {batch}")
  dacs <- h5read(hdf5_file, name = paste0("/Batches/", batch, "/Dacs")) %>% 
    data.table() %>% 
    setnames("dacs")
  
  # Read Ref_to_signal (information of number of dacs mapping to each ref position)
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
  
  log_success("Finished processing {batch}")
  return(Ref_to_signal)
}


###############################################################################
# Process signal mappings
###############################################################################

# Append reference mappings
add_mapping_to_dacs <- function(dacs, mappings, type){
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
    ][
      , dacs_norm := (dac - mean(dac))/sd(dac), by = read_id
    ]
}

# Calculate mean difference and u-test value
calculate_current_diff <- function(dt1, dt2, min_u_val, min_cov = 3){
  dt1[
      dt2, on = "contig_id", nat := i.nat
    ][
      , mean_nat := pblapply(
          nat,
          function(x){mean(x$dacs_norm, na_rm = TRUE)},
          cl = ifelse(parallel, arg$threads, NULL)) %>% 
        unlist()
    ][
      , mean_pcr := pblapply(
          pcr, 
          function(x){mean(x$dacs_norm, na_rm = TRUE)},
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
        FUN = function(nat, pcr, n_nat_map, n_pcr_map){
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
        cl = arg$threads
      )
    ][
      , u_val_weighted := rolling_mean(u_val, n = arg$u_val_weight_window, weigth_dropoff = arg$u_val_weight_dropoff)
    ][
      , u_val_weighted_log := rolling_mean_log(u_val, n = arg$u_val_weight_window, weigth_dropoff = arg$u_val_weight_dropoff)
    ][
      , event := u_val < min_u_val
    ][
      , event_weighted := u_val_weighted_log < min_u_val*1e3
    ]
}

###############################################################################
# Visualisation of events in one contig
###############################################################################
# Plots 

plot_events <- function(dt){
  guides_format <- guides(
      fill = guide_legend(title = "Event"),
      size = guide_legend(title = "Event"),
      col = guide_legend(title = "Event"),
      alpha = "none"
    )
  plot_col <- c("#99712d", "#2d5599")
  geom_point_default <- geom_jitter(alpha = 0.1, size = 0.1)
  
  
  p_1 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    ggplot(aes(x = contig_index, y = mean_dif)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = mean_dif), y = 0, size = 0.2) +
    geom_point(aes(fill = event_weighted,  size = event_weighted), shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (Mean difference)"
    )  + 
    scale_size_manual(values = c(0.4, 2)) +
    scale_fill_manual(values = plot_col) +
    guides_format
  
  p_2 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    ggplot(aes(x = contig_index, y = u_val_weighted_log)) +
    geom_hline(yintercept = 1) +
    geom_segment(aes(x = contig_index, xend = contig_index, yend = u_val_weighted_log), y = 0, size = 0.2) +
    geom_point(aes(fill = event_weighted, size = event_weighted), shape = 21) +
    labs(
      x = "Contig Position",
      y = "NAT vs. PCR (U-value)"
    ) +
    scale_size_manual(values = c(0.4, 2)) +
    scale_y_log10() +
    scale_fill_manual(values = plot_col) +
    guides_format 
    
  
  p_3 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    ggplot(aes(x = n_pcr, y = n_nat, col = event_weighted, fill = event_weighted)) +
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
    guides_format
  
  p_4 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    ggplot(aes(x = n_pcr_map, y = n_nat_map, col = event, fill = event_weighted)) +
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
    guides_format
  
  p_5 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    mutate(
      mappings = n_pcr_map / n_nat_map
    ) %>% 
    ggplot(aes(y = mappings, x = u_val_weighted_log, col = event_weighted)) +
    labs(
      x = "U-value",
      y = "nr. PCR / nr. NAT"
    ) + 
    geom_point_default +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format
  
  p_6 <- dt %>% 
    filter(!is.na(event_weighted)) %>% 
    ggplot(aes(x = mean_dif, y = u_val_weighted_log, col = event_weighted)) +
    labs(
      x = "Mean difference",
      y = "U-value"
    ) + 
    geom_point_default +
    scale_y_log10() +
    scale_color_manual(values = plot_col) +
    guides_format
  
  p_all <- (p_1 | p_2) / (p_3 | p_4) / (p_5 | p_6) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Significant events in ", signal_mappings$contig %>% unique() %>% `[`(1)," (U-value threshold: ", arg$min_u_val, ")")
    )  &
    theme_bw() 
  return(p_all)
}

# Importing reads to reference mappings
mapping_pcr <- fread("/shared-nfs/SH/samples/zymo/megalodon/pcr_test/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4, V2) %>% 
  setnames(c("read_id", "contig", "start", "direction")) %>% 
  mutate(
    direction = ifelse(direction == 16, "fwd", "rev")
  )

mapping_nat <- fread("/shared-nfs/SH/samples/zymo/megalodon/nat_test2/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4, V2) %>% 
  setnames(c("read_id", "contig", "start", "direction")) %>% 
  mutate(
    direction = ifelse(direction == 16, "fwd", "rev")
  )

###############################################################################
# PCR reads
###############################################################################

if (!file.access(glue("{arg$out}/singal_mappings_pcr.tsv"), mode = 4) == 0 | arg$overwrite) {
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
  
  signal_mappings_pcr <- add_mapping_to_dacs(dacs_pcr, mapping_pcr, "pcr") 
  
  # Save signal mappings
  fwrite(
    signal_mappings_pcr,
    glue("{arg$out}/singal_mappings_pcr.tsv")
  )
  
  signal_mappings_pcr <- signal_mappings_pcr %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "pcr")
    

} else {
  signal_mappings_pcr <- fread(glue("{arg$out}/singal_mappings_pcr.tsv")) %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "pcr")
}

###############################################################################
# NAT reads
###############################################################################

if (!file.access(glue("{arg$out}/singal_mappings_nat.tsv"), mode = 4) == 0 | arg$overwrite) {
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
          reads_ids = mapping_nat[contig %in% arg$contig_plot,][, read_id][1:3415]
        )
      }
    ) %>% 
    rbindlist()
  
  h5closeAll()
  rm(hdf5_nat)
  # Add mapping mapping to reference
  signal_mappings_nat <- add_mapping_to_dacs(dacs_nat, mapping_nat, "nat") 
  
  # Save signal mappings
  fwrite(
    signal_mappings_nat,
    glue("{arg$out}/singal_mappings_nat.tsv")
  )
  
  signal_mappings_nat <- signal_mappings_nat %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "nat")
  
} else {
  signal_mappings_nat <- fread(glue("{arg$out}/singal_mappings_nat.tsv")) %>% 
    group_nest_dt(contig_id, contig_index, contig, direction, .key = "nat")
}

###############################################################################
# Calculate current difference
###############################################################################

if (!file.access(glue("{arg$out}/current_difference.tsv"), mode = 4) == 0 | arg$overwrite) {
  current_difference <- calculate_current_diff(
      signal_mappings_pcr, 
      signal_mappings_nat, 
      arg$min_u_val
    )
  current_difference[, `:=`(nat = NULL, pcr = NULL)]
   
  # Save processed signal mappings
  fwrite(
    current_difference,
    glue("{arg$out}/current_difference.tsv")
  )
} else {
  current_difference <- fread(glue("{arg$out}/current_difference.tsv"))
}

p_all <- plot_events(current_difference)

ggsave(
  width = unit(8, "cm"),
  height = unit(10, "cm"),
  p_all, 
  filename = glue("{arg$out}/current_difference_{arg$contig_plot}.png"), 
  device = "png"
)

###############################################################################
# Evaluation of PCR batch 1 vs batch 0
###############################################################################

signal_mappings_pcr_0 <- add_mapping_to_dacs(dacs_pcr[batch == "Batch_0", ], mapping_pcr, "pcr") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "pcr")

signal_mappings_pcr_1 <- add_mapping_to_dacs(dacs_pcr[batch == "Batch_1", ], mapping_pcr, "nat") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "nat")

calculate_current_diff(signal_mappings_pcr_0, signal_mappings_pcr_1, arg$min_u_val)


plot_events(signal_mappings_pcr_0[
  n_nat_map > 5 & n_pcr_map > 5
])

ggsave(
  width = unit(8, "cm"),
  height = unit(10, "cm"),
  p_pcr_vs_pcr, 
  filename = glue("{arg$out}/current_difference_pcr_vs_pcr_{arg$contig_plot}.png"), 
  device = "png"
)

###############################################################################
# Evaluation of differenct U-value cut offs
###############################################################################

# pcr vs. pcr and nat vs. pcr plot function
plot_PvP_NvP <- function(thresholds, col){
  thresholds %>% 
  lapply(
    function(x){
      n_NvP <- sum(current_difference[, ..col] < x, na.rm = TRUE)
      n_PvP <- sum(current_difference_pcr[, ..col] < x, na.rm = TRUE)
      list(
        x,
        n_NvP,
        n_NvP / nrow(na.omit(current_difference[, ..col])), 
        n_PvP, 
        n_PvP / nrow(na.omit(current_difference_pcr[, ..col]))
      )
    }) %>% 
  rbindlist() %>% 
  setnames(paste0("V", 1:5), c("U_val", "n_NvP", "p_NvP", "n_PvP", "p_PvP")) %>%  
  tidyr::pivot_longer(starts_with("p_"), values_to = "proportion", names_to = "type") %>% 
  tidyr::pivot_longer(starts_with("n_"), values_to = "count", names_to = "type2") %>% 
  filter(substr(type, 2, 4) == substr(type2, 2, 4)) %>%  
  ggplot() +
  aes(x = U_val, y = proportion, fill = type) +
  geom_line() +
  geom_point(shape = 21, size = 1.5) +
  scale_x_log10() +
  scale_y_log10(labels = percent) +
  labs(
    x = "P-value threshold",
    y = "Proportion of events"
  ) +
  theme_bw()
}

# p-value
p2.1 <-  10^-seq(-1, 20, by = 0.5) %>% 
  plot_PvP_NvP(col = "u_val") +
  labs(title = "P-value")

# p-value weighted
p2.2 <-  10^-seq(1, 20, by = 0.5) %>% 
  plot_PvP_NvP(col = "u_val_weighted") +
  labs(title = glue("P-value weighted, dropoff = {arg$dropoff}"))

# p-value log weighted 
p2.3 <-  10^-seq(1, 20, by = 0.5) %>% 
  plot_PvP_NvP(col = "u_val_weighted_log") +
  labs(title = glue("P-value log weighted, dropoff = {arg$dropoff}"))
  
# Mean difference cutoff
p2.4 <- 10^-seq(-1, 8, by = 0.5) %>% 
  plot_PvP_NvP("mean_dif") +
  xlab("Mean difference threshold") +
  labs(title = "Mean difference")


p2.1 + p2.2 /
p2.3 + p2.4
plot_layout(guides = "collect") +
  plot_annotation(
    title =  "Estimation of event cut-off"
  ) & 
  theme_bw()



plot_events(current_difference[contig_index %in% 7.15e4:7.25e4, ][
  n_nat_map > 5 & n_pcr_map > 5
])



