log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = snakemake@params[["data_dir"]],
  molecule.coordinates = FALSE
)

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")
snakemake@source("../../scripts/_standard_seurat_analysis/_add_metadata_post_seurat_load_xenium.R")

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
