#!/usr/bin/env Rscript

.libPaths(c("/shared-nfs/SH/code/Rlib", .libPaths()))

setwd(here::here())
pacman::p_load(
  "ape",
  "data.table",
  "argparse"
)
# Read in arguments
parser <- ArgumentParser()

parser$add_argument("-e", "--embedding",  nargs=1,
                    help="File containing low dimentional representation of methylation features")
parser$add_argument("-a", "--assembly",  nargs=1,
                    help="Assembly fasta file")
parser$add_argument("-o", "--out",  nargs=1,
                    help="Output folder [Default: next to profile path]")


args <- parser$parse_args()

if (is.null(args$embedding)) {
  stop(sprintf("Required argument \"-e, --embedding\" is missing"))
}
if (is.null(args$assembly)) {
  stop(sprintf("Required argument(\"-a, --assembly\") is missing"))
}
# Check assembly input
if (file.access(args$assembly) == -1) {
  stop(sprintf("Assembly file does not exist (\"%s\")", args$assembly))
}
# Check embedding input
if (file.access(args$embedding) == -1) {
  stop(sprintf("Profile file does not exist (\"%s\")", args$embedding))
}
# Set out to default of none specified
if (is.null(args$out)) {
  args$out <- dirname(args$embedding)
}


source("/shared-nfs/SH/code/scripts/functions.R")

# args <-list()
# args$assembly <- "/shared-nfs/SH/samples/zymo/assembly_100x.fasta"
# args$embedding <- "/shared-nfs/SH/samples/zymo/methbin/AE/binnning_AE_representation.tsv"
# args$out <- "/shared-nfs/SH/samples/zymo/methbin/AE"
features <- fread(args$embedding)

# Fasta sequences of contigs
fasta <- read.fasta(args$assembly)

# Low dimentional representation of contigs
embedding <- fread(args$embedding) %>%
  setnames(paste0("dim", 1:length(1:length(.)))) %>%
  cbind(contig = features$contig) %>%
  cbind(setnames(as.data.table(Rtsne(.)$Y), c("tsne1", "tsne2")))  %>%
  left_join(get_GC_percent(fasta), by = "contig") %>%
  mutate(
    cluster = hdbscan(subset(.,select =  grepl("dim", names(.))), minPts = 2)$cluster
  )

# Plotting
ggplot(embedding, aes(x=tsne1, y=tsne2, size = contig_length, fill = as.character(cluster))) +
  geom_point(shape = 21) +
  scale_size(range = c(0.1, 10))

# Write clusters
write_cluster_fasta(
  df = embedding,
  metagenome_fasta = fasta,
  path_out = paste0(args$out, "/bins_AE/")
)
