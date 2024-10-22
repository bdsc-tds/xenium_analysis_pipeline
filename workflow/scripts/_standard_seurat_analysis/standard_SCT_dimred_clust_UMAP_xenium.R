# not sure we need to store log-normalized and scaled data for the Xenium assay. It will take quire some space and im not sure its used in any of the downstream analyses
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)

default_assay  <- snakemake@params[["default_assay"]]
n_dims         <- snakemake@params[["n_dims"]]
resolution     <- snakemake@params[["resolution"]]

dims           <- 1:n_dims

# Read post QC Seurat
xe             <- readRDS(snakemake@input[["xe"]])

# Run standard SCTransform and downstream analysis
xe             <- xe %>% SCTransform(assay = default_assay) %>% RunPCA(npcs = n_dims) %>% FindNeighbors(dims = dims) %>% FindClusters(resolution = resolution)
xe             <- xe %>% RunUMAP(dims = dims)

# `seurat_clusters` gets owerwritten every time `Seurat::FindClusters()` is called. Save it in a new "proterted" column
xe@meta.data   <- xe@meta.data %>% mutate(SCT_seurat_clusters = seurat_clusters)

# Save parameters in seurat object
xe@misc$standard_seurat_analysis_meta <- list(
  n_dims     = n_dims,
  dims       = dims,
  resolution = resolution
)

# Save post QC seurat
saveRDS(
  xe, 
  file = file.path(snakemake@output[["xe"]])
)
