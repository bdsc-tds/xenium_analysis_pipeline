log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)
library(puRCTD)

snakemake@source("../../scripts/utils/readwrite.R")


xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])
purified_counts <- Read10X_h5(snakemake@input[["fully_purified_counts"]])
purified_counts_metadata <- read_parquet(snakemake@input[["fully_purified_counts_metadata"]])

###### Balancing purification level ######

# Create purified Seurat
xe_purified <- CreateSeuratObject(counts = purified_counts, assay = "Xenium", meta.data = purified_counts_metadata)
rm("purified_counts")

### Run `spot_class_balanced` purification (ie, purify doublets and do not touch singlets)
### Done by balancing raw and purified datasets
xe <- AddMetaData(xe, rctd@results$results_df)
xe_balanced_spot_class <- balance_raw_and_purified_data_by_spot_class(
  xe_raw = xe,
  xe_purified = xe_purified,
  assay = "Xenium", # should be param, but can wait 
  spot_class_key = "spot_class",
  DO_swap_lables = TRUE  # should be param, but can wait 
)

# Output spot-class-based purified counts
write10xCounts(
  path = snakemake@output[["corrected_counts"]], 
  x = GetAssayData(xe_balanced_spot_class, assay = "Xenium", layer = "counts")
)

write_parquet(
  xe_balanced_spot_class@meta.data %>% select(all_of(colnames(purified_counts_metadata))), 
  snakemake@output[["corrected_counts_metadata"]]
)
