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

################################## PURIFICATION ############################################
message("Runing SPLIT::split...\n")
### Run full purification (ie., purify all cells except for highly confident singlets (ie., those that do not have the second cell type)
res_split <- SPLIT::split(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'),
  rctd = rctd,
  DO_purify_singlets = TRUE
)
message("Done SPLIT::split...\n")

message("Saving SPLIT::split output ...\n")
# Output purified counts
write10xCounts(path = snakemake@output[["corrected_counts"]], x = res_split$purified_counts)
write_parquet(as.data.frame(res_split$cell_meta), snakemake@output[["corrected_counts_metadata"]])
message("Done saving SPLIT::split output ...\n")