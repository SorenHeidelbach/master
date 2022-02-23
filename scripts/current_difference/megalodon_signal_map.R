#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/software/Rlib" , .libPaths()))
setwd(here::here())
pacman::p_load(
  "data.table",
  "dplyr",
  "ggplot2", 
  "rhdf5",
  "seqinr",
  "microbenchmark",
  "foreach",
  "Rsamtools",
  "stringr",
  "patchwork",
  "logger"
)



group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}
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

###############################################################################
# Importing signal_mapping.HDF5 
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

  dacs_test <- 
    mapply(
      function(len, end, id, ref, pos){
        list(
          dac <- dacs$dacs[(end-len+1):end],
          dac_pos = (end-len+1):(end),
          read_id = rep(id, len),
          ref = rep(ref, len),
          pos = rep(pos, len)
        )
      },
      Ref_to_signal$n_dac,
      Ref_to_signal$start,
      Ref_to_signal$read_id,
      Ref_to_signal$ref,
      Ref_to_signal$pos,
      
      SIMPLIFY = FALSE
    ) %>% 
    rbindlist()
  log_success("Finished processing {batch}")
  return(dacs_test)
}

mapping_pcr <- fread("/shared-nfs/SH/samples/zymo/megalodon/pcr_test/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

mapping_nat <- fread("/shared-nfs/SH/samples/zymo/megalodon/nat_test2/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

log_info("Processing PCR mappings")
hdf5_pcr <- H5Fopen(arg$signal_mapping_pcr)
dacs_pcr <- h5ls(hdf5_pcr) %>% filter(group == "/Batches") %>% pull(name) %>% 
  mclapply(
    mc.cores = 2,
    function(batch){
      load_mapping_hdf5(
        hdf5_pcr, 
        batch = batch,
        reads_ids = mapping_pcr[contig %in% "lf_contig3",][, read_id]
      )
    }
  ) %>% 
  rbindlist()

log_info("Processing NAT mappings")
hdf5_nat <- H5Fopen(arg$signal_mapping_nat)
dacs_nat <- h5ls(hdf5_nat) %>% filter(group == "/Batches") %>% pull(name) %>% 
  mclapply(
    mc.cores = 2,
    function(batch){
      load_mapping_hdf5(
        hdf5_nat, 
        batch = batch,
        reads_ids = mapping_nat[contig %in% "lf_contig3",][, read_id]
  )
    }
  ) %>% 
  rbindlist()

# Plot of normalised dacs of a read
dacs_nat[read_id == unique(dacs_nat$read_id)[1],][
    , dacs_norm := (V1 - mean(V1))/sd(V1), by = read_id
  ] %>% 
  ggplot(aes(x = pos, y = dacs_norm)) +
  geom_line() +
  theme_minimal()

# Load in reference mappings
add_mapping_to_dacs <- function(dacs, mappings, type){
  dacs[
      mappings, on = "read_id", `:=`(contig = i.contig, start = i.start)
    ][
      , type := type
    ][
      , pos := pos - 1
    ][
      , contig_index := start + pos
    ][
      , contig_id := paste0(contig, "_", contig_index)
    ][
  , dacs_norm := (V1 - mean(V1))/sd(V1), by = read_id
]
}
signal_mappings <- add_mapping_to_dacs(dacs_pcr, mapping_pcr, "pcr") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "pcr")

signal_mappings2 <- add_mapping_to_dacs(dacs_nat, mapping_nat, "nat") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "nat")

arg$min_u_val <- 1e-8
signal_mappings[
    signal_mappings2, on = "contig_id", nat := i.nat
  ][
    , mean_nat := lapply(nat, function(x){mean(x$dacs_norm, na_rm = TRUE)}) %>% unlist()
  ][
    , mean_pcr := lapply(pcr, function(x){mean(x$dacs_norm, na_rm = TRUE)}) %>% unlist()
  ][
    , mean_dif := mean_nat - mean_pcr
  ][
    , n_nat := lapply(nat, nrow) %>% unlist()
  ][
    , n_pcr := lapply(pcr, nrow) %>% unlist()
  ][
    , u_val := mapply(FUN = function(nat, pcr){
      if(length(nat$dacs_norm) > 0 & length(pcr$dacs_norm) > 0) {
        wilcox.test(nat$dacs_norm, pcr$dacs_norm)$p.value
      } else {
        NA
      }},
      nat,
      pcr
    )
    ][
      , u_val_weighted := rolling_mean(u_val, n = arg$u_val_weight_window, weigth_dropoff = arg$u_val_weight_dropoff)
    ][
      , u_val_weighted_log := rolling_mean_log(u_val, n = arg$u_val_weight_window, weigth_dropoff = arg$u_val_weight_dropoff)
    ][
      , event := u_val < min_u_val
    ][
      , event_weighted := u_val_weighted_log < min_u_val*1e3
    ][
      , n_nat_map := lapply(nat, function(x) x$read_id %>% unique %>% length) %>% unlist()
    ][
      , n_pcr_map := lapply(pcr, function(x) x$read_id %>% unique %>% length) %>% unlist()
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
    scale_size_manual(values = c(0.9, 3)) +
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
      x = "PCR current values at contig position",
      y = "NAT current values at contig position"
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
      x = "Number of PCR reads mapped to contig position",
      y = "Number of NAT reads mapped to contig position"
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
      y = "Proportion of PCR to NAT reads mapped"
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
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

mapping_nat <- fread("/shared-nfs/SH/samples/zymo/megalodon/nat_test2/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

# Loading signal mappings
log_info("Processing PCR mappings")
hdf5_pcr <- H5Fopen(arg$signal_mapping_pcr)
dacs_pcr <- h5ls(hdf5_pcr) %>% filter(group == "/Batches") %>% pull(name) %>% 
  lapply(
    function(batch){
      load_mapping_hdf5(
        hdf5_pcr, 
        batch = batch,
        reads_ids = mapping_pcr[contig %in% "lf_contig3",][, read_id]
      )
    }
  ) %>% 
  rbindlist()

log_info("Processing NAT mappings")
hdf5_nat <- H5Fopen(arg$signal_mapping_nat)
dacs_nat <- h5ls(hdf5_nat) %>% filter(group == "/Batches") %>% pull(name) %>% 
  lapply(
    function(batch){
      load_mapping_hdf5(
        hdf5_nat, 
        batch = batch,
        reads_ids = mapping_nat[contig %in% "lf_contig3",][, read_id]
      )
    }
  ) %>% 
  rbindlist()


signal_mappings <- add_mapping_to_dacs(dacs_pcr, mapping_pcr, "pcr") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "pcr")

signal_mappings2 <- add_mapping_to_dacs(dacs_nat, mapping_nat, "nat") %>% 
  group_nest_dt(contig_id, contig_index, contig, .key = "nat")

calculate_current_diff(signal_mappings, signal_mappings2, arg$min_u_val)

# Save processed signal mappings
fwrite(
  signal_mappings[
    , unlist(nat, recursive = FALSE), 
    by = eval(names(signal_mappings)[!(names(signal_mappings) %in% c("pcr", "nat"))])],
  glue("{arg$out}/current_difference.tsv")
)

p_all <- plot_events(signal_mappings)

ggsave(
  width = unit(8, "cm"),
  height = unit(10, "cm"),
  p_all, 
  filename = glue("{arg$out}/current_difference_{arg$contig_plot}.png"), 
  device = "png"
)






