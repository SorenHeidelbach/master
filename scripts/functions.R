
# Setting costum library path
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
if (!require("pacman")) install.packages("pacman", lib = "/shared-nfs/SH/code/Rlib")
if (!require("BiocManager")) install.packages("BiocManager", lib = "/shared-nfs/SH/code/Rlib")

# Loading libraries with pacman
pacman::p_load(
  "tidyr",
  "dplyr",
  "seqinr"
)
##-----------------------------------------------------------------------------
get_GC_percent <-  function(fasta, min_size = 0) {
  # Calculate GC% from fasta file path
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

normalise_cols <- function(numeric_dataframe){
  # Normalise values in columns of data frame
  numeric_dataframe %>%
    apply(MARGIN = 2, function(x) (x-min(x)) / (max(x) - min(x))) %>%
    as.data.frame()
}

standardise_cols <- function(numeric_dataframe){
  # Standardise values in columns of data frame
  numeric_dataframe %>%
    apply(MARGIN = 2, function(x) (x - mean(x))/sd(x)) %>%
    as.data.frame()
}

group_nest_dt <- function(dt, ..., .key = "data"){
  # goup_nest version for data.table
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

unnest_dt <- function(dt, col, id){
  # unnest version for data.table 
  stopifnot(is.data.table(dt))
  by <- substitute(id)
  col <- substitute(unlist(col, recursive = FALSE))
  dt[, eval(col), by = eval(by)]
}

set_sum_1 <- function(x){
  x/sum(x)
}

generate_pascal_trianble <- function(levels = 1000){
  # Generate list with all layers of pascal triangle to specified level (e.g. 1, 11, 121, 1331, 14641)
  x <- 1
  pascal <- list()
  for (i in 1:levels) {x <- c(0, x) + c(x, 0); pascal[[i]] <- x }
  return(pascal)
}

embed_vector <- function(vec, pascal){
  # Turn size n vector into size 5 vector using pascal triangle
  if(length(vec) > 1000) vec <- vec[sample(1:length(vec), 1000) %>% sort]
  if(length(vec) <=5 ) return(vec)
  level <- length(vec) - 5
  vec_out <- lapply(
    1:5,
    function(i){
      sum(vec[i:(level+i)] * pascal[[level]]) / (2^level)
    }
  ) %>%  unlist()
  return(vec_out)
}

rolling_mean <- function(x, n = 2, weight_dropoff = 0.75){
  # Calculate rolling mean with distance weigth decay
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

mean_vec <- function(vec, frequency = 5){
  if (length(vec) %% frequency != 0) {
    warning("One of the vectors is of not correct length, remainder values are dropped")
    vec <- vec[1:(length(vec)-(length(vec) %% frequency))]
  }
  
  matrix <- matrix(vec, ncol = frequency, byrow = TRUE)
  mean_vec <- (matrix %>% apply(2, sum)) / nrow(matrix)
  return(mean_vec)
}
mean_vec_dist <- function(vec1, vec2){
  # Calculate mean vector and get distance between these
  if (length(vec1) %% freq != 0 | length(vec2) %% freq != 0) {
    warning("One of the vectors is of not correct length, remainder values are dropped")
    vec1 <- vec1[1:length(vec1)-(length(vec1) %% freq)]
    vec2 <- vec2[1:length(vec2)-(length(vec2) %% freq)]
  }
  matrix1 <- matrix(vec1, ncol = freq, byrow = TRUE)
  matrix2 <- matrix(vec2, ncol = freq, byrow = TRUE)
  data.frame(
    dist = rbind(
      (matrix1 %>% apply(2, sum)) / nrow(matrix1),
      (matrix2 %>% apply(2, sum)) / nrow(matrix2)
    ) %>% dist() %>% c()
  )
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
