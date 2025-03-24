log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)
library(SPLIT)

snakemake@source("../../scripts/utils/readwrite.R")


xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])
sp_neigh_df <- read_parquet(snakemake@input[["spatial_neighborhood_scores"]])
purified_counts <- Read10X_h5(snakemake@input[["fully_purified_counts"]])
purified_counts_metadata <- read_parquet(snakemake@input[["fully_purified_counts_metadata"]])

###### Balancing purification level ######

# Create purified Seurat
xe_purified <- CreateSeuratObject(counts = purified_counts, meta.data = purified_counts_metadata)
rm("purified_counts")

### Run `score_based` purification (ie, purify doublets singlets that have the spatial neighborhood score larges than the threshold)
### Done by balancing raw and purified datasets
xe <- AddMetaData(xe, rctd@results$results_df)
xe <- AddMetaData(xe, sp_neigh_df)

xe_balanced_score <- balance_raw_and_purified_data_by_score(
  xe_raw = xe,
  xe_purified = xe_purified,
  assay = "Xenium", # should be param, but can wait 
  threshold = snakemake@params[["score_threshold"]],
  score_name = snakemake@params[["score_name"]],
  DO_swap_lables = TRUE # should be param, but can wait 
)

# Output score-based purified counts
write10xCounts(
  path = snakemake@output[["corrected_counts"]], 
  x = GetAssayData(xe_balanced_score, assay = "Xenium", layer = "counts")
)

write_parquet(
  xe_balanced_score@meta.data %>% select(all_of(colnames(purified_counts_metadata))), 
  snakemake@output[["corrected_counts_metadata"]]
)
