log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Export the Xenium query of RCTD for rctd-py: counts as 10x v3 HDF5, cell coordinates as parquet.
# The query is the same object and assay rctd.R reads.

library(Seurat)
library(arrow)

snakemake@source("../../../scripts/utils/readwrite.R")

query_path <- snakemake@input[["query"]]
xe_assay   <- snakemake@params[["xe_assay"]]

# Load the query and check the assay
xe <- readRDS(query_path)
if (!xe_assay %in% Assays(xe)) stop("Assay '", xe_assay, "' is not found in the xenium object.")

# Counts layer, as in rctd.R
counts <- GetAssayData(xe, assay = xe_assay, layer = "counts")
message(sprintf("Query: %d genes x %d cells", nrow(counts), ncol(counts)))

dir.create(dirname(snakemake@output[["counts"]]), recursive = TRUE, showWarnings = FALSE)
if (file.exists(snakemake@output[["counts"]])) file.remove(snakemake@output[["counts"]])
write10xCounts(
  path    = snakemake@output[["counts"]],
  x       = counts,
  type    = "HDF5",
  genome  = "unknown",
  version = "3"
)

# Cell ids and coordinates (ST_1 / ST_2, the coordinates rctd.R gives SpatialRNA)
metadata <- data.frame(
  cell_id = colnames(xe),
  ST_1    = xe@meta.data[colnames(xe), "ST_1"],
  ST_2    = xe@meta.data[colnames(xe), "ST_2"]
)
write_parquet(metadata, snakemake@output[["metadata"]])
