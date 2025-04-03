log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)
library(SPLIT)

snakemake@source("../../scripts/utils/readwrite.R")

message("Reading data...")

xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])
tr_neigh_df <- read_parquet(snakemake@input[["transcriptomic_neighborhood_scores"]])
purified_counts <- Read10X_h5(snakemake@input[["fully_purified_counts"]])

md <- open_dataset(snakemake@input[["fully_purified_counts_metadata"]], thrift_string_size_limit = 1000000000) # force reading of large files
purified_counts_metadata <- md %>%
  collect() %>%
  mutate(across(everything(), unname))

message("Done reading data.")

###### Balancing purification level ######
message("Creating seurat objects...")
# Create purified Seurat
xe_purified <- CreateSeuratObject(counts = purified_counts, assay = "Xenium", meta.data = purified_counts_metadata)
rm("purified_counts")

### Run `spot_class_balanced` purification (ie, purify doublets and do not touch singlets)
### Done by balancing raw and purified datasets
xe <- AddMetaData(xe, rctd@results$results_df)
xe <- AddMetaData(xe, tr_neigh_df)
message("Done creating seurat objects.")

message("Balancing by spot class")

xe_balanced_spot_class <- balance_raw_and_purified_data_by_spot_class(
  xe_raw = xe,
  xe_purified = xe_purified,
  default_assay = "Xenium", # should be param, but can wait 
  spot_class_key = "spot_class",
  DO_swap_lables = TRUE  # should be param, but can wait 
)
message("Done balancing by spot_class")

message("Writing data...")

# Output spot-class-based purified counts
write10xCounts(
  path = snakemake@output[["corrected_counts"]], 
  x = GetAssayData(xe_balanced_spot_class, assay = "Xenium", layer = "counts")
)

write_parquet(
  xe_balanced_spot_class@meta.data %>% select(all_of(c(colnames(purified_counts_metadata), "swap"))), 
  snakemake@output[["corrected_counts_metadata"]]
)
message("Done writing data...")