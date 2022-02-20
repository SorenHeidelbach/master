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
  H5Fclose(hdf5_file)
  
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

hdf5_pcr <- H5Fopen(arg$signal_mapping_pcr)
test <- h5ls(hdf5_pcr) %>% filter(group == "/Batches") %>% pull(name) %>% 
  lapply(
    function(batch){
      load_mapping_hdf5(
        hdf5_pcr, 
        batch = batch
      )
    }
  )
dacs_pcr <- load_mapping_hdf5(
  arg$signal_mapping_pcr, 
  batch = arg$batch,
  reads_ids = mapping_pcr %>% subset(contig == "bs_contig1") %>% pull(read_id) %>% unique()
  )

dacs_nat <- load_mapping_hdf5(arg$signal_mapping_nat, reads_ids = mapping_nat %>% subset(contig == "bs_contig1") %>% pull(read_id) %>% unique())

# Plot of normalised dacs of a read
dacs_nat[read_id == unique(dacs_nat$read_id)[1],][
    , dacs_norm := (V1 - mean(V1))/sd(V1), by = read_id
  ] %>% 
  ggplot(aes(x = pos, y = V1)) +
  geom_line()

# Load in reference mappings
signal_mappings <- rbind(
  dacs_pcr %>% 
    left_join(mapping_pcr, by = "read_id") %>% 
    mutate(
      pos = pos - 1,
      contig_index = start + pos,
      contig_id = paste0(contig, "_", contig_index),
      type = "pcr"
    ),
  dacs_nat %>% 
    left_join(mapping_nat, by = "read_id") %>% 
    mutate(
      pos = pos - 1,
      contig_index = start + pos,
      contig_id = paste0(contig, "_", contig_index),
      type = "nat"
    )
) 
signal_mappings[
  , dacs_norm := (V1 - mean(V1))/sd(V1), by = read_id
]

signal_mappings <- signal_mappings %>%  group_nest_dt(contig_id, contig_index, contig) 

# p_dac_norm <- signal_mappings[contig_id %in% paste0("bs_contig1_4990", 10:20), ] %>% 
#   ggplot(aes(x = as.character(contig_index), y = dacs_norm, fill = type)) +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(shape = 21, size = 0.5, position = position_jitterdodge(seed = 2)) +
#   theme_bw() +
#   theme(axis.text.x.bottom = element_text(angle = 90))
# 
# p_dac <- signal_mappings[contig_id %in% paste0("bs_contig1_4990", 10:20), ] %>% 
#   ggplot(aes(x = as.character(contig_index), y = V1, fill = type)) +
#   geom_boxplot(outlier.alpha = 0) +
#   geom_jitter(shape = 21, size = 0.5, position = position_jitterdodge(seed = 2), ) +
#   theme_bw() +
#   theme(axis.text.x.bottom = element_text(angle = 90))
# 
# p_dac_norm / p_dac
library(doParallel)  # will load parallel, foreach, and iterators
cl <- makeCluster(30)
registerDoParallel(cl)
signal_mappings2 <- signal_mappings %>% 
  left_join(signal_mappings %>%
    mutate(
      chunk = floor(contig_index/1000)
    ) %>% 
    split(f = .$chunk) %>% 
    mclapply( 
      mc.cores = 30,
      function(x){
        apply(
          x,
          MARGIN = 1,
          function(dt){
            pcr = dt$data[type == "pcr",]
            nat = dt$data[type == "nat",]
            if ((nrow(pcr) != 0) & (nrow(nat) != 0)){
              p <- wilcox.test(pcr$dacs_norm, nat$dacs_norm)$p.value
              list(
                contig_id = dt$contig_id,
                p_val = p,
                dacs_nat_minus_pcr = mean(nat$dacs_norm) - mean(pcr$dacs_norm)
              )
            }else {
              list(
                contig_id = dt$contig_id,
                p_val = NA,
                dacs_nat_minus_pcr = mean(nat$dacs_norm) - mean(pcr$dacs_norm)
              )
            }
          }
        ) %>% 
          rbindlist()
      }
    ) %>% 
      rbindlist()
)

arg$min_u_val <- 1e-8

p_mean_dif_event <- signal_mappings2 %>% 
  filter(!is.na(p_val)) %>% 
  mutate(
    event = ifelse(p_val < arg$min_u_val, TRUE, FALSE)
  ) %>% 
  ggplot(aes(x = contig_index, y = dacs_nat_minus_pcr)) +
  geom_hline(yintercept = 0) +
  geom_segment(aes(x = contig_index, xend = contig_index, yend = dacs_nat_minus_pcr), y = 0, size = 0.2) +
  geom_point(aes(fill = event, size = event), shape = 21) +
  scale_size_manual(values = c(0.9, 3)) +
  guides(
    fill = guide_legend(title = "Event"),
    size = guide_legend(title = "Event")
  ) +
  labs(
    x = "Contig Position",
    y = "NAT vs. PCR (Mean difference)"
  ) #+ ggforce::facet_zoom(xlim = c(218500, 219500))

p_pval_event <- signal_mappings2 %>% 
  filter(!is.na(p_val)) %>% 
  mutate(
    event = ifelse(p_val < arg$min_u_val, TRUE, FALSE)
  ) %>% 
  ggplot(aes(x = contig_index, y = p_val)) +
  geom_hline(yintercept = 1) +
  geom_segment(aes(x = contig_index, xend = contig_index, yend = p_val), y = 0, size = 0.2) +
  geom_point(aes(fill = event, size = event), shape = 21) +
  scale_size_manual(values = c(0.9, 3)) +
  scale_y_log10() +
  guides(
    fill = guide_legend(title = "Event"),
    size = guide_legend(title = "Event")
  ) +
  labs(
    x = "Contig Position",
    y = "NAT vs. PCR (U-value)"
  ) #+ggforce::facet_zoom(xlim = c(218500, 219500))
p_work <- (p_pval_event | p_mean_dif_event) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Significant events in ", contig_select, " (U-value threshold: ", arg$min_u_val, ")")
  )  &
  theme_bw()
p_work




