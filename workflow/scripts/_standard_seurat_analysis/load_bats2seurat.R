log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

snakemake@source("../../scripts/utils/run_time_utils.R")

library(arrow)
library(dplyr)

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = snakemake@params[["data_dir"]],
  molecule.coordinates = FALSE
)

expected_counts <- data.frame(read_parquet(snakemake@input[["expected_counts"]]))
rownames(expected_counts) <- expected_counts$cell_id
expected_counts <- expected_counts %>% select(-cell_id)
expected_counts <- expected_counts[, colnames(expected_counts)[!grepl(
  snakemake@params[["control_gene_pat"]],
  colnames(expected_counts)
)]]

xe <- replace_counts_in_seurat(
  xe,
  expected_counts
)

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")
snakemake@source("../../scripts/_standard_seurat_analysis/_add_metadata_post_seurat_load_xenium.R")

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
