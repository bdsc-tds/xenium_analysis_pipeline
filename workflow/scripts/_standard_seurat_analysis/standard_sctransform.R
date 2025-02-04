log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)
library(arrow)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

xe <- readRDS(snakemake@input[[1]])

xe <- xe %>% SCTransform(assay = snakemake@params[["default_assay"]])

xe@misc$standard_seurat_analysis_meta <- list(
  normalisation_id = snakemake@params[["normalisation_id"]]
)

normalised_data <- GetAssayData(
  xe,
  assay = snakemake@params[["normalised_assay"]],
  layer = snakemake@params[["normalised_layer"]]
)

# Save post QC seurat
saveRDS(
  xe,
  file = file.path(snakemake@output[["obj"]])
)

write_parquet(
  data.frame(
    cell = colnames(normalised_data)
  ),
  sink = file.path(snakemake@output[["cells"]])
)

write_parquet(
  as.data.frame(t(normalised_data)),
  sink = file.path(snakemake@output[["counts"]])
)
