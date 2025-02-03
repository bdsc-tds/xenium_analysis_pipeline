log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Annotate Xenium from Chromium using SingleR

library(Seurat)
library(dplyr)
library(SingleR)
library(arrow)

# Load common reference-based parameters
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_header.R")


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'data' 
test_layer <- 'counts' 

# Load reference and query data
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_load_data.R") # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

## Make sure chrom data are log-normalized
chrom               <- NormalizeData(chrom)

CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]

# Generate reference object
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_generate_reference_obj.R")

ref_labels <- chrom@meta.data %>% pull(annotation_level)

# SingleR parameters (ideally, should be provided with `_other_options` from the experiment.yml or config.yml)
selected_genes <- snakemake@params[["genes"]] %||% "de"
de_method <- snakemake@params[["de_method"]] %||% "t"
de_n <- snakemake@params[["de_n"]] %||% 25
aggr_ref <- snakemake@params[["aggr_ref"]] %||% TRUE
aggr_args <- snakemake@params[["aggr_args"]] %||% list(rank = 50, power = 0.7)

# Annotate with SingleR
message("Running SingleR...")
singler_result <- SingleR(
  test = GetAssayData(xe, assay = xe_assay, layer = test_layer),
  ref = GetAssayData(chrom, assay = ref_assay, layer = ref_layer),
  labels = ref_labels,
  genes = selected_genes,
  de.method = de_method,
  de.n = de_n,
  aggr.ref = aggr_ref,
  aggr.args = aggr_args
)

singler_scores                <- singler_result$scores %>% as.data.frame()
singler_scores$cell_id        <- colnames(xe)
singler_scores                <- singler_scores %>% select(cell_id, everything())

singler_labels                <- data.frame(cell_id = colnames(xe), label = singler_result$labels)

# Save annotation
saveRDS(singler_result, snakemake@output[[1]])
write_parquet(singler_labels, snakemake@output[[2]])
write_parquet(singler_scores, snakemake@output[[3]])


