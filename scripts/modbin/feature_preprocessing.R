#!/usr/bin/env Rscript
.libPaths(c("/shared-nfs/SH/code/Rlib", .libPaths()))

setwd(here::here())
pacman::p_load(
  "ape",
  "data.table",
  "argparse",
  "tictoc"
)
# Read in arguments
parser <- ArgumentParser()

parser$add_argument("-p", "--profile",  nargs=1,
                    help="Path to Nanodisco methylation profile")
parser$add_argument("-a", "--assembly",  nargs=1,
                    help="Assembly fasta file")
parser$add_argument("-o", "--out",  nargs=1,
                    help="Output folder [Default: next to profile path]")


parser$add_argument("--min_contig_size", default=10000, type="integer",
                    help="Minimum of size of contig to considered for binning [Default: %(default)s]")
parser$add_argument("--filter_method", nargs=1, default="max",
                    help="Filtration method for removal of low info features [Default: %(default)s]")
parser$add_argument("--missing_value_imputation", nargs=1, default=0,
                    help="Method for filling missing values [Default: %(default)s]")
parser$add_argument("--min_current_distance", default=1.5, type="double",
                    help="Minimum value for current distance score [Default: %(default)s]")
parser$add_argument("--min_num_occurance", default=5, type="double",
                    help="Minimum number of mapped reads tyo current distance position [Default: %(default)s]")

args <- parser$parse_args()

if (is.null(args$profile)) {
  stop(sprintf("Required argument \"-p, --profile\" is missing"))
}
if (is.null(args$assembly)) {
  stop(sprintf("Required argument(\"-a, --assembly\") is missing"))
}
# Check assembly input
if (file.access(args$assembly) == -1) {
  stop(sprintf("Assembly file does not exist (\"%s\")", args$assembly))
}
# Check profile input
if (file.access(args$profile) == -1) {
  stop(sprintf("Profile file does not exist (\"%s\")", args$profile))
}
# Set out to default of none specified
if (is.null(args$out)) {
  args$out <- dirname(args$profile)
}


source("/shared-nfs/SH/code/scripts/functions.R")
# Load methylation features
##-----------------------------------------------------------------------------
timer=c(Sys.time(), Sys.time())
# Load profile 
message("Loading profile")
profile <- readRDS(args$profile)
profile <- as.data.table(profile)
difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()

message("  Metainformation")
# Contig meta information
contigs <- profile %>%
  attributes() %>%
  `[[`("contig_coverage") %>%
  filter(contig_length > args$min_contig_size)  %>%
  arrange(-contig_length) %>%
  mutate(size_order = as.numeric(row.names(.))) %>% 
  mutate(contig=chr) %>% 
  select(-chr)
difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()

message("  Removing small contigs")
# Combine contig information with profile data frame
profile <- profile[
    contigs, on = "contig", (names(contigs)[-7]) := mget(paste0("i.", names(contigs)[-7]))
  ][
    contig_length >= args$min_contig_size
  ]
difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()
##-----------------------------------------------------------------------------
# Profile into features

message("Filtering")
# Removing feature vectors that doesn't contain any informative values
if(args$filter_method == "max"){
  message("  Keeping all values in feature vector if just one of the values are above minimum")
  profile <- profile[
    , c("max_nb_occurrence", "max_dist_score") := list(max(nb_occurrence), max(dist_score)), by = .(motif, distance_motif)
  ][
    max_nb_occurrence > args$min_num_occurance & max_dist_score > args$min_current_distance, 
  ]
  difftime(Sys.time(), timer[2])
  timer[2] <- Sys.time()
} else {
  message("  Choose valid filtering method")
  stop()
}

message("Indentifying missing feature values")
# pviot longer to get each feature in a column
profile <- profile[
    , feature_name := paste0(motif, "_", distance_motif)
  ][
    , .(contig, feature_name, dist_score)
  ][
    expand.grid(contig = unique(contig), feature_name = unique(feature_name)), on = .(contig, feature_name)
  ]

difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()

# Setting missing values according to specified method
contigs <- profile$contig
if (args$missing_value_imputation %in% "mean") {
  message("  Mising values are filled with feature vector mean")
  profile[
    , dist_score := ifelse(is.na(dist_score), mean(dist_score, na.rm = TRUE), dist_score), by = feature_name
  ]
} else if (is.numeric(args$missing_value_imputation)) {
  message(paste0("  Missing values are replaced with ", args$missing_value_imputation))
  profile[
    , dist_score := ifelse(is.na(dist_score), args$missing_value_imputation, dist_score), by = feature_name
  ]
} else if (args$missing_value_imputation == "NA") {
  message("  Missing values are kept as NA")
  # returns dataframe with NA values
} else {
  message("  Please set a valid 'missing_value' ('mean' or replacement number or 'NA')")
  stop()
}


profile <- profile %>%
  dcast(contig ~ feature_name, value.var = "dist_score")

message("  Missing values treated")
difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()


message("  Writing features to tsv")
dir.create(args$out, showWarnings = FALSE)
fwrite(profile, file = paste0(args$out, "/features.tsv"), sep="\t")

difftime(Sys.time(), timer[2])
timer[2] <- Sys.time()
