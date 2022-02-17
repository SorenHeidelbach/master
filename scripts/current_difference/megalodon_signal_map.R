#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
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
  "stringr"
)

  

arg <- list()
arg$reference <- "/shared-nfs/SH/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta"
arg$signal_mapping_pcr <- "/shared-nfs/SH/samples/zymo/megalodon/pcr_test/signal_mappings.hdf5"
arg$signal_mapping_nat <- "/shared-nfs/SH/samples/zymo/megalodon/nat_test2/signal_mappings.hdf5"

arg$batch <- "Batch_0"

###############################################################################
# Importing signal_mapping.HDF5 
###############################################################################

# Load in reference
ref <- read.fasta(arg$reference)

load_mapping_hdf5 <- function(arg_mapping, reads_ids = NA){
  # Open HDF5
  hdf5_file <- H5Fopen(arg_mapping)
  
  # h5ls(hdf5_file)
  
  # Read dacs (current measurements)
  dacs <- h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Dacs")) %>% 
    data.table() %>% 
    setnames("dacs")
  
  # Read Ref_to_signal (information of number of dacs mapping to each ref position)
  Ref_to_signal <- h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Ref_to_signal")) %>% 
    data.table() %>% 
    setnames("ref_to_signal")
  
  # Read Reference (Interger representation of read bases)
  Reference = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Reference")) %>% 
    data.table() %>% 
    setnames("reference")
  
  # Read read information (mainly ID, dacs length, Ref/sig length)
  reads <- h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Dacs_lengths")) %>% 
    data.table() %>% 
    setnames("Dacs_lengths") %>% 
    cbind("Ref_to_signal_lengths" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Ref_to_signal_lengths"))) %>% 
    cbind("Reference_lengths" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/Reference_lengths"))) %>% 
    #cbind("digitisation" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/digitisation"))) %>% 
    cbind("offset" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/offset"))) %>% 
    cbind("range" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/range"))) %>% 
    #cbind("scale_frompA" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/scale_frompA"))) %>% 
    #cbind("shift_frompA" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/shift_frompA"))) %>% 
    apply(MARGIN = 2, FUN = as.numeric) %>% 
    data.table(stringsAsFactors = FALSE) %>% 
    cbind("read_id" = h5read(hdf5_file, name = paste0("/Batches/", arg$batch, "/read_id")))
  H5Fclose(hdf5_file)
  
  print("asd2")
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
  
  
  print("asd3")
  if (!is.na(reads_ids)) {
    Ref_to_signal <- Ref_to_signal %>% 
      subset(read_id %in% reads_ids)
  }
  tictoc::tic()
  dacs_test <- 
    mcmapply(
      mc.cores = 1,
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
  return(dacs_test)
}

# Sanity check for similarity between raw file and dacs
p_current_v_dac <- list()
p_current_v_dac <- lapply(
  reads$read_id,
  function(x){
    # Dacs from signal mapping
    print(x)
    dac_plot <- dacs_test[
        read_id == x, 
      ][
        , pos := 1:.N
      ]
    list(
      ggplot() +
        geom_line(
          data = dac_plot,
          mapping = aes(x = pos, y = V1, group = 1, ), 
          size = 0.5
          ) +
        # geom_line(
        #   data = raw_plot,
        #   mapping = aes(x = pos, y = signal+1000, group = 1), 
        #   size = 0.5
        # ) +
        theme_bw() +
        ggforce::facet_zoom(x = c(1.8e4, 2e4)),
      ggplot(Ref_to_signal[read_id == x, ], aes(x = pos, y = signal_frame)) +
        geom_line( 
          size = 0.5
        ) +
        theme_bw() 
    )
  }
)

group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

mapping_pcr <- fread("/shared-nfs/SH/samples/zymo/megalodon/pcr_test/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

mapping_nat <- fread("/shared-nfs/SH/samples/zymo/megalodon/nat_test2/mappings_sorted_view.txt") %>% 
  select(V1, V3, V4) %>% 
  setnames(c("read_id", "contig", "start"))

dacs_pcr <- load_mapping_hdf5(arg$signal_mapping_pcr, reads_ids = mapping_pcr %>% subset(contig == "bs_contig1") %>% pull(read_id) %>% unique())
dacs_nat <- load_mapping_hdf5(arg$signal_mapping_nat, reads_ids = mapping_nat %>% subset(contig == "bs_contig1") %>% pull(read_id) %>% unique())


# Load in reference mappings

signal_mappings <- rbind(
  dacs_pcr %>% 
    left_join(mapping, by = "read_id") %>% 
    mutate(
      pos = pos - 1,
      contig_index = paste0(contig, "_", start + pos),
      type = "pcr"
    ),
  dacs_nat %>% 
    left_join(mapping, by = "read_id") %>% 
    mutate(
      pos = pos - 1,
      contig_index = paste0(contig, "_", start + pos),
      type = "nat"
    )
) %>% group_nest_dt(contig_index)

utest <- signal_mappings$data %>% 
  lapply(
    function(dt){
      pcr = dt[type == "pcr",]
      nat = dt[type == "nat",]
      if ((nrow(pcr) != 0) & (nrow(nat) != 0)){
        wilcox.test(nat$V1, pcr$V1)
      }
      
    }
  )


