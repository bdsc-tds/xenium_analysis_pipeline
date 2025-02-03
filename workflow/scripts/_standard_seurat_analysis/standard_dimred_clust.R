log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)
library(tibble)
library(arrow)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

default_assay <- snakemake@params[["default_assay"]]
n_dims <- snakemake@params[["n_dims"]]
resolution <- snakemake@params[["resolution"]]

dims <- 1:n_dims

# Read post QC Seurat
xe <- readRDS(snakemake@input[[1]])

xe <- xe %>% RunPCA(npcs = n_dims) %>%
    FindNeighbors(dims = dims) %>%
    FindClusters(resolution = resolution) %>%
    RunUMAP(dims = dims)

# `seurat_clusters` gets overwritten every time `Seurat::FindClusters()` is called. Save it in a new "protected" column
xe@meta.data <- xe@meta.data %>% mutate(SCT_seurat_clusters = seurat_clusters)

# Save parameters in seurat object
xe@misc$standard_seurat_analysis_meta <- c(
  xe@misc$standard_seurat_analysis_meta,
  list(
    n_dims = n_dims,
    dims = dims,
    resolution = resolution
  )
)

# Save results
saveRDS(
  xe, 
  file = file.path(snakemake@output[["obj"]])
)

write_parquet(
  as.data.frame(
    Embeddings(xe, reduction = "pca")
  ) %>% rownames_to_column("cell_id"),
  sink = file.path(snakemake@output[["pca"]])
)

write_parquet(
  as.data.frame(
    Embeddings(xe, reduction = "umap")
  ) %>% rownames_to_column("cell_id"),
  sink = file.path(snakemake@output[["umap"]])
)
