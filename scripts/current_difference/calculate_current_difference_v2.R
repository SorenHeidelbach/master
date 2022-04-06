#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/software/Rlib" , .libPaths()))
pacman::p_load(
  "data.table",
  "dplyr",
  "rhdf5",
  "seqinr",
  "stringr",
  "logger",
  "glue",
  "pbapply",
  "parallel",
  "R.utils",
  "argparse",
  "here"
)
setwd(here())
source("scripts/functions.R")
print("ONLY E. coli plasmid is evaluated currently!!!! hardcoded in ")
###############################################################################
# Utility functions
###############################################################################


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
  unlink(paste0(arg$out, "/temp"))
}

###############################################################################
# Input arguments
###############################################################################
# Read in arguments
parser <- ArgumentParser()

#parser$add_argument("--data_dir", help='Directionary containing output from "calculate_current_difference.R"', required = TRUE)

parser$add_argument("--signal_mapping_nat",  nargs=1, help="Signal mappings for natural reads", required = TRUE)
parser$add_argument("--signal_mapping_pcr",  nargs=1, help="Signal mappings for PCR reads", required = TRUE)
parser$add_argument("--read_mapping_nat",  nargs=1, help="Read mappings for natural reads (bam | sort | view)", required = TRUE)
parser$add_argument("--read_mapping_pcr",  nargs=1, help="Read mappings for PCR reads (bam | sort | view)", required = TRUE)
parser$add_argument("--overwrite", action="store_true", help = "Overwrite previous results in out folder")
parser$add_argument("--pcr_vs_pcr", action="store_true", help = "Calculate PCR vs. PCR")
parser$add_argument("--out",  nargs=1, help="Output folder [default data_dir]")
parser$add_argument("-t", "--threads", default = 10, type = "integer", help = "Threads to use")
parser$add_argument("--min_mappings", default = 3, type = "integer", help = "Minimum coverage of position to be included (both NAT and PCR)")

arg <- parser$parse_args()
arg$min_p_val <- 1e-9

arg <- list()
arg$signal_mapping_pcr <- "/shared-nfs/SH/data/zymoHMW/PCR/megalodon/megalodon/pcr_HMW-zymo_Z2_continued/signal_mappings.hdf5"
arg$signal_mapping_nat <- "/shared-nfs/SH/data/zymoHMW/NAT/megalodon/signal_mappings.hdf5"
arg$read_mapping_pcr <- "/shared-nfs/SH/data/zymoHMW/PCR/megalodon/megalodon/pcr_HMW-zymo_Z2_continued/mappings.view.txt.gz"
arg$read_mapping_nat <- "/shared-nfs/SH/data/zymoHMW/NAT/megalodon/mappings.view.txt.gz"
arg$out <- "/shared-nfs/SH/results/zymoHMW/current_difference/test1_20220404"
arg$contig_plot <- c("bs_contig1", "lf_contig1")
arg$overwrite <- FALSE
arg$save_intermediate_files <- TRUE
arg$threads <- 20
arg$min_mappings <- 3

# arg$out <- "/shared-nfs/SH/results/zymoHMW/current_difference/test1_20220331"
arg$p_val_weight_window <- 3 # val*2 + 1
arg$p_val_weight_dropoff <- 0.9


dir.create(arg$out, showWarnings = FALSE)
jsonlite::write_json(arg, path = glue::glue("{arg$out}/args.json"), pretty = TRUE)



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
    cbind("scale_frompA" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/scale_frompA"))) %>% 
    cbind("shift_frompA" = h5read(hdf5_file, name = paste0("/Batches/", batch, "/shift_frompA"))) %>% 
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
    , unlist(dac, recursive = FALSE), by = eval(names(Ref_to_signal)[!(names(Ref_to_signal) %in% c("dac"))])
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
  min_p_val, 
  min_cov = 3, 
  data_splits = 1,
  parallel = TRUE){
  # log_info("Setting up cluster")
  # cl <- makeForkCluster(nnodes = arg$threads)
  pboptions(nout = data_splits)
  log_info("Distance and Mean difference")
  dt1[
      dt2, on = .(contig_id, direction), nat := i.nat
    ][
      , dist := dist(rbind(nat, pcr)) %>% c(), by = .(contig_id, direction)
    ][
      , mean_dif := mean(nat) - mean(pcr), by = .(contig_id, direction)
    ]
  log_info("Wilcox test")
  dt2 <- dt1[
      , .(p_val = wilcox.test(nat, pcr)$p.value), by = .(contig_id, dist, mean_dif, contig, contig_index, direction)
    ]
  log_info("Weighting p-values")
  dt2[
      , p_val_weighted_log := rolling_mean_log(p_val, contig_index = contig_index, n = arg$p_val_weight_window, weight_dropoff = arg$p_val_weight_dropoff), by = contig
    ][
      , p_val_bonf := p.adjust(p_val, n = length(p_val), method = "bonferroni")
    ][
      , p_val_bonf_weighted_log := rolling_mean_log(p_val_bonf, contig_index = contig_index, n = arg$p_val_weight_window, weight_dropoff = arg$p_val_weight_dropoff), by = contig
    ][
      , p_val_BH := p.adjust(p_val, n = length(p_val), method = "BH")
    ][
      , p_val_BH_weighted_log := rolling_mean_log(p_val_BH, contig_index = contig_index, n = arg$p_val_weight_window, weight_dropoff = arg$p_val_weight_dropoff), by = contig
    ]
}

###############################################################################
##                                                                           ##
##                                                                           ##
##                     Processing                                            ##
##                                                                           ##
##                                                                           ##
###############################################################################

log_info("Loading read mappings")
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

if (!file.access(glue("{arg$out}/singal_mappings_pcr.tsv.gz"), mode = 4) == 0 | arg$overwrite) {
  # Loading signal mappings
  log_info("Processing PCR signal mappings")
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
          reads_ids = mapping_pcr$read_id[mapping_pcr$contig == "Escherichia_coli_plasmid"]
        )
      }
    ) %>% 
    rbindlist()
  h5closeAll()
  rm(hdf5_pcr)
  
  
  pascal <- generate_pascal_trianble()
  signal_mappings_pcr <- signal_mappings_pcr_unnested[
        n_dac >= 1, .(dac_embed = as.numeric(embed_vector(dac, pascal))), by = .(read_id, contig_id, contig_index, contig, direction)
      ][
        , .(pcr = mean_vec(dac_embed)), by = .(contig_id, contig_index, contig, direction)
      ]
  
  # Save signal mappings
  log_info("Saving processed PCR signal mappings")
  fwrite(signal_mappings_pcr_unnested, glue("{arg$out}/singal_mappings_pcr.tsv.gz"), compress = "gzip")
  
  log_success("Succesfully processed PCR signal mappings")
} else {
  log_info("Loading PCR signal mappings")
  signal_mappings_pcr_unnested <- fread(glue("{arg$out}/singal_mappings_pcr.tsv.gz"))
}

signal_mappings_pcr <- signal_mappings_pcr_unnested %>% 
  group_nest_dt(contig_id, contig_index, contig, direction, .key = "pcr")
rm(signal_mappings_pcr_unnested)

###############################################################################
# NAT reads
###############################################################################

if (!file.access(glue("{arg$out}/singal_mappings_nat.tsv.gz"), mode = 4) == 0 | arg$overwrite) {
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
          reads_ids = mapping_nat$read_id[mapping_nat$contig == "Escherichia_coli_plasmid"]
        )
      }
    ) %>% 
    rbindlist()
  h5closeAll()
  rm(hdf5_nat)
  
  
  pascal <- generate_pascal_trianble()
  signal_mappings_nat <- signal_mappings_nat_unnested[
      n_dac >= 1, .(dac_embed = as.numeric(embed_vector(dac, pascal))), by = .(read_id, contig_id, contig_index, contig, direction)
    ][
      , .(nat = mean_vec(dac_embed)), by = .(contig_id, contig_index, contig, direction)
    ]
  
  
  log_info("Saving processed NAT signal mappings")
  fwrite(signal_mappings_nat_unnested, glue("{arg$out}/singal_mappings_nat.tsv.gz"), compress = "gzip")
  
  log_success("Succesfully processed NAT signal mappings")
} else {
  log_info("Loading NAT signal mappings")
  signal_mappings_nat_unnested <- fread(glue("{arg$out}/singal_mappings_nat.tsv.gz"))
}

signal_mappings_nat <- signal_mappings_nat_unnested2 %>% 
  group_nest_dt(contig_id, contig_index, contig, direction, .key = "nat")

signal_mappings_nat <- signal_mappings_nat_unnested %>% 
  group_nest_dt(contig_id, contig_index, contig, direction, .key = "nat")
rm(signal_mappings_nat_unnested)

###############################################################################
# Calculate current difference
###############################################################################

if (!file.access(glue("{arg$out}/current_difference.tsv.gz"), mode = 4) == 0 | arg$overwrite) {
  log_info("Calculating current difference and statistics")
  current_difference <- copy(signal_mappings_pcr)
  current_difference <- calculate_current_diff(
    current_difference, 
    signal_mappings_nat, 
    arg$min_p_val
  )
  current_difference[, `:=`(nat = NULL, pcr = NULL)]
  setorderv(current_difference, c("contig", "contig_index"))
  
  # Save processed signal mappings
  log_info("Saving current difference")
  fwrite(
    current_difference,
    compress = "gzip"
  )
  log_success("Performed currrent difference")
}

rm(signal_mappings_nat)
rm(signal_mappings_pcr)

###############################################################################
# Evaluation of PCR vs. PCR
###############################################################################
if ((!file.access(glue("{arg$out}/current_difference_pcr.tsv.gz"), mode = 4) == 0 | arg$overwrite) & arg$pcr_vs_pcr ) {
  signal_mappings_pcr_unnested <- fread(glue("{arg$out}/singal_mappings_pcr.tsv.gz"))
  
  # Split sample in two by contig/ref
  pcr_read_id_1 <- mapping_pcr[mapping_pcr[, .I[ read_id %in% sample(read_id, length(read_id)/2) ], by=contig]$V1]$read_id
  pcr_read_id_2 <- mapping_pcr[! read_id %in% pcr_read_id_1, ]$read_id
  
  signal_mappings_pcr_0 <- signal_mappings_pcr_unnested[read_id %in%  pcr_read_id_1, ][
      n_dac >= 1, .(dac_embed = as.numeric(embed_vector(dac, pascal))), by = .(read_id, contig_id, contig_index, contig, direction)
    ][
      , .(pcr = mean_vec(dac_embed)), by = .(contig_id, contig_index, contig, direction)
    ]
  signal_mappings_pcr_1 <- signal_mappings_pcr_unnested[read_id %in% pcr_read_id_2, ][
      n_dac >= 1, .(dac_embed = as.numeric(embed_vector(dac, pascal))), by = .(read_id, contig_id, contig_index, contig, direction)
    ][
      , .(nat = mean_vec(dac_embed)), by = .(contig_id, contig_index, contig, direction)
    ]
  
  # Current difference
  current_difference_pcr <- copy(signal_mappings_pcr_0)
  current_difference_pcr <- calculate_current_diff(current_difference_pcr, signal_mappings_pcr_1,  min_p_val = arg$min_p_val)
  
  # Save
  current_difference_pcr[, `:=`(nat = NULL, pcr = NULL)]
  
  fwrite(
    current_difference_pcr,
    glue("{arg$out}/current_difference_pcr.tsv.gz"),
    compress = "gzip"
  )
}


