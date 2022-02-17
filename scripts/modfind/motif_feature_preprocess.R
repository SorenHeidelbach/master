#!/usr/bin/Rscript
.libPaths( c( "/shared-nfs/SH/code/Rlib" , .libPaths()))
setwd(here::here())
pacman::p_load(
  "data.table",
  "dplyr",
  "argparse",
  "jsonlite",
  "ggplot2"
)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("-d", "--dif",  nargs=1,
                    help="Path to current difference")
parser$add_argument("-o", "--out",  nargs=1,
                    help="Output folder")

parser$add_argument("--event_frame_size", default=16, type="double",
                    help="Size of frame included in features after event detection")
parser$add_argument("--event_u_val_threshold", default=1e-5, type="double",
                    help="Threshold value for U-test p-value")
args <- parser$parse_args()

if (is.null(args$dif)) {
  stop(sprintf("Required argument \"-d, --dif\" is missing"))
}
if (is.null(args$out)) {
  stop(sprintf("Required argument(\"-o, --out\") is missing"))
}

write_json(args, paste0(args$out, "/preprocess_settings.json"), pretty=TRUE)
# Variables
dif_path <- args$dif
dif_feature_out <- paste0(args$out, "/event_features.tsv") 
dif_metainfo_out <- paste0(args$out, "/metainfo.tsv")

# Hyperparameters
event_frame <-args$event_frame_size
event_frame <- c(-ceiling(event_frame/2):-1, 1:ceiling(event_frame/2)) # Rounded up if uneven 

event_threshold <- args$event_u_val_threshold

# Read in current difference of metagenome
print("Reading in difference RDS file")
dif <- data.table(readRDS(dif_path))

# Assign 1 to failed U-tests
print("Assigning value to failed U-tests")
dif <- dif[
    is.na(u_test_pval), u_test_pval := 1
  ][
    , event := 0
  ][
    u_test_pval <= event_threshold, event := 1
  ]
print("Expanding with adjacent current values within event frame")
# Add adjacent current difference values to events
dif <- dif %>% setorder(dir, contig, position)
dif <- dif[
    , c("strand", "t_test_pval", "N_wga", "N_nat") := rep(list(NULL), 4)
  ][
    , paste0("rel_pos_", event_frame) := lapply(event_frame, function(x) shift(mean_diff, -x)), by = list(contig, dir)
  ]

# Add corresponding reverse current difference values
print("Joining forwards and reverse values")
dif <- merge(
  dif[dir=="fwd"], dif[dir=="rev"], 
  all = TRUE, 
  by = c("contig", "position"),
  suffixes = c("fwd", "rev"))

# Remove non event rows

dif_events <- dif$eventrev |  dif$eventfwd
dif_features <- dif[dif_events,]

# Plot contigs with high event overlap 
# dif_features %>% 
#   group_by(contig) %>% 
#   mutate(
#     dist_pos = sort(position) - shift(sort(position), type = "lead"),
#     dist_pre = shift(sort(position)) - sort(position),
#     overlap = dist_pos < max(event_frame) | dist_pre < max(event_frame)
#   ) %>% 
#   summarise(proportion_overlap = sum(overlap)/nrow(.)) %>% 
#   filter(proportion_overlap %in% tail(sort(proportion_overlap), 20)) %>% 
#   mutate(
#     contig = factor(contig, levels = contig[order(-proportion_overlap)])
#   ) %>% 
#   ggplot(aes(y = proportion_overlap, x=contig)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x.bottom = element_text(angle = 90))


# Save meta information
print("Saving feature metadata")
fwrite(dif_features[, c("contig", "position")], dif_metainfo_out, sep = "\t")

print("Reformatting data.table")
dif_features <- dif_features %>% 
  select(-contig, -position, -starts_with("dir"), -starts_with("u_test_pval"), -starts_with("event")) %>% 
  as.data.frame() %>% 
  # Make sure sequence order is preserved
  `[`(, 
      c(cbind(
          # Postions of forward reads
          c(
            2:(ceiling(ncol(.)*1/4)), # upstream                     
            1, # even tposition
            ((ceiling(ncol(.)*1/4)+1)):(ceiling(ncol(.)*2/4)) # downstream
            ),
          # Positions of reverse reads
          c(((ceiling(ncol(.)*2/4))+2):(ceiling(ncol(.)*3/4)), 
            ( ceiling(ncol(.)*2/4))+1, 
            ( ceiling(ncol(.)*3/4) +1):ncol(.)
            )
        ))
  )

# Rowwise normalisation
print("Normalisation")
dif_features <- dif_features %>% 
  apply(
    MARGIN = 2,
    function(x) {
      return((x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    }) %>% 
  rbind() %>% 
  as.data.frame() 
# Random value pertubation
dif_features[is.na(dif_features)] <- runif(sum(is.na(dif_features)), min = 0, max= 0.05)

# Save features
print("Saving features")
fwrite(dif_features, dif_feature_out, sep = "\t")

# Plot visualsing event in contig region
# dif %>% 
#   filter(contig == "contig_1" & position > 1.7*10e4 & position < 1.8*10e4) %>% 
#   mutate(
#     mean_difffwd = ifelse(is.na(mean_difffwd), 0, mean_difffwd),
#     eventfwd = ifelse(eventfwd==1, "True", "False")
#     ) %>% 
#   ggplot(aes(x = position, y = mean_difffwd, fill = as.character(eventfwd))) +
#   geom_errorbar(aes(ymax = mean_difffwd, x = position), ymin = 0, alpha = 0.6) +
#   geom_point(aes(size = eventfwd), shape = 21) +
#   scale_size_manual(values = c(1, 3), guide = 'none') +
#   guides(fill = guide_legend(title = "Event")) +
#   scale_fill_manual(values = c("skyblue3", "tomato1")) +
#   labs(x = "Contig Position", y = "Mean Current  Difference") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0.05, 0.05)) +
#   theme_bw() 
