# Annotate Xenium from Chromium using SingleR

library(Seurat)
library(dplyr)
library(SingleR)

# Load common reference-based parameters
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_header.R")) 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'data' 
test_layer <- 'counts' 

# Load reference and query data
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_load_data.R")) # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

# Generate reference object
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_generate_reference_obj.R")) 

ref_labels <- chrom@meta.data %>% pull(annotation_level)

# SingleR parameters (ideally, should be provided with `_other_options` from the experiment.yml or config.yml)
selected_genes <- snakemake@params[["singleR_genes"]] %||% "de"
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

singler_scores                <- singler_result$scores
singler_labels                <- singler_result$labels
rownames(singler_scores)      <- colnames(xe)
names(singler_labels)         <- colnames(xe)

# Save annotation
saveRDS(singler_result, snakemake@output[[1]])
write.csv(singler_labels, snakemake@output[[2]])
write.csv(singler_scores, snakemake@output[[3]])


