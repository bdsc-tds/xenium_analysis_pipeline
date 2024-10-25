log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)

convert2numeric <- function(val) {
  if (is.numeric(val)) return(val)

  ret <- as.numeric(val)

  stopifnot(!is.na(ret))

  return(ret)
}

default_assay <- snakemake@params[["default_assay"]]
default_layer <- snakemake@params[["default_layer"]]

# Read raw seurat
xe <- readRDS(snakemake@input[[1]])

# QC thresholds, either global or gene panel specific
min_counts <- convert2numeric(snakemake@params[["min_counts"]])
min_features <- convert2numeric(snakemake@params[["min_features"]])
max_counts <- convert2numeric(snakemake@params[["max_counts"]])
max_features <- convert2numeric(snakemake@params[["max_features"]])
min_cells <- convert2numeric(snakemake@params[["min_cells"]])

# Dims of raw data 
n_cells_raw <- ncol(xe)
n_genes_raw <- nrow(xe)

# Rename metadata to simplify access for filtering
xe@meta.data[["nFeature"]] <- xe@meta.data[[paste0("nFeature", "_", default_assay)]]
xe@meta.data[["nCount"]] <- xe@meta.data[[paste0("nCount", "_", default_assay)]]

# Filter Cells
xe <- xe %>% subset(
  (nFeature > min_features & nFeature < max_features) &
    (nCount > min_counts & nCount < max_counts)
)
message(paste("N =", n_cells_raw - ncol(xe), " cells were removed"))

# Filter Genes 
gene_cells <- rowSums(GetAssayData(xe, assay = default_assay, layer = default_layer) > 0)
genes_to_keep <- gene_cells[gene_cells > min_cells] %>% names()

xe <- xe[genes_to_keep,]
message(paste("N =", n_genes_raw - nrow(xe), " genes were removed"))

# Add QC thresholds to seurat object to keep track 
xe@misc$QC_metadata <- list(
  min_counts = min_counts,
  min_features = min_features,
  max_counts = max_counts,
  max_features = max_features,
  min_cells = min_cells
)

# Save post QC seurat
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
