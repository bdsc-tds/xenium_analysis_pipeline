log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(arrow)

snakemake@source("../../scripts/utils/run_time_utils.R")

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- readRDS(snakemake@input[["raw_obj"]])
corrected_counts <- Read10X_h5(snakemake@input[["corrected_counts"]])

# Substitute original reads with corrected ones

original_cell_ids  <- colnames(xe)
corrected_cell_ids <- colnames(corrected_counts)

if(sum(!corrected_cell_ids %in% original_cell_ids) > 0){
  warning("Cell ids in corrected count have been modofied, reqeuire `cell_id` column in metadata to map ids")
}

xe <- replace_counts_in_seurat(
  xe = xe, 
  new_counts = t(corrected_counts), 
  cell_coords = xe@meta.data[corrected_cell_ids, c("ST_1", "ST_2")]
)

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")
snakemake@source("../../scripts/_standard_seurat_analysis/_add_metadata_post_seurat_replace_counts.R")

# Save object
saveRDS(
  xe,
  file = file.path(snakemake@output[[1]])
)
