log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)
library(arrow)

xe <- readRDS(snakemake@input[[1]])

xe <- xe %>% NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData()

xe@misc$standard_seurat_analysis_meta <- list(
  normalisation_id = snakemake@params[["normalisation_id"]]
)

data <- GetAssayData(
  xe,
  assay = snakemake@params[["assay"]],
  layer = snakemake@params[["data_layer"]]
)

scale.data <- GetAssayData(
  xe,
  assay = snakemake@params[["assay"]],
  layer = snakemake@params[["scale_data_layer"]]
)

# Save post QC seurat
saveRDS(
  xe,
  file = file.path(snakemake@output[[1]])
)

write_parquet(
  data.frame(
    cell = colnames(data)
  ),
  sink = file.path(snakemake@output[["cells"]])
)

write_parquet(
  as.data.frame(t(data)),
  sink = file.path(snakemake@output[["data"]])
)

write_parquet(
  as.data.frame(t(scale.data)),
  sink = file.path(snakemake@output[["scale_data"]])
)
