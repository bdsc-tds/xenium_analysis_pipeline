log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

default_assay <- snakemake@params[["default_assay"]]

xe <- readRDS(snakemake@input[[1]])

xe <- xe %>% SCTransform(assay = default_assay)

xe@misc$standard_seurat_analysis_meta <- list(
  normalisation_id = snakemake@params[["normalisation_id"]]
)

# Save post QC seurat
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
