log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Annotate Xenium from Chromium using Seurat Transfer

library(Seurat)
library(dplyr)
library(arrow)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

# Load common reference-based parameters
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_header.R")


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'data' 
test_layer <- 'data' 

# Load reference and query data
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_load_data.R") # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

## Make sure chrom data are log-normalized
DefaultAssay(chrom) <- ref_assay
# chrom               <- NormalizeData(chrom)


DefaultAssay(xe)    <- xe_assay
#xe                  <- NormalizeData(xe)

CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]

# Generate reference object
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_generate_reference_obj.R")

# ? Should chrom data be re-normalized besed on shared genes???

ref_labels <- chrom@meta.data %>% pull(annotation_level) 

xe_chrom_common_genes <- intersect(rownames(xe), rownames(chrom))

normalization.method <- xe@misc$standard_seurat_analysis_meta$normalisation_id
if (normalization.method == "sctransform") {
    normalization.method <- "SCT"
} else if (normalization.method == "lognorm") {
    normalization.method <- "LogNormalize"
} else {
    stop("Unknown normalization method: ", normalization.method)
}

# Annotate with Seurat transfer 

dims <- snakemake@params[["min_dim"]]:snakemake@params[["max_dim"]]

anchors       <- FindTransferAnchors(
    reference = chrom,
    query = xe,
    normalization.method = normalization.method,
    features = xe_chrom_common_genes,
    dims = dims
)
predictions   <- TransferData(anchorset = anchors, refdata = ref_labels, dims = dims)
xe            <- AddMetaData(xe, metadata = predictions)


scores                <- predictions %>% select(!c("predicted.id", "prediction.score.max"))
colnames(scores)      <- gsub("prediction.score.", "", colnames(scores))
labels                <- predictions$predicted.id

# Convert result to data.frames for parquet output
scores$cell_id        <- colnames(xe)
scores                <- scores %>% select(cell_id, everything())

labels_df             <- data.frame(cell_id = colnames(xe), label = labels)

# Save annotation
saveRDS(predictions, snakemake@output[[1]])
write_parquet(labels_df, snakemake@output[[2]])
write_parquet(scores, snakemake@output[[3]]) 

