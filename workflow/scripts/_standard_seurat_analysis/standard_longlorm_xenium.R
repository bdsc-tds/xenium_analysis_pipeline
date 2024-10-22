# not sure we need to store log-normalized and scaled data for the Xenium assay. It will take quire some space and im not sure its used in any of the downstream analyses
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)

default_assay  <- snakemake@params[["default_assay"]]

xe             <- readRDS(snakemake@input[["xe"]])

xe             <- xe %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

# Save post QC seurat
saveRDS(
  xe, 
  file = file.path(snakemake@output[["xe"]])
)
