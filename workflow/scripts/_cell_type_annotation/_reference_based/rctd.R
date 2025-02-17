log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Annotate Xenium from Chromium using RCTD

library(Seurat)
library(dplyr)
library(spacexr)
library(arrow)
library(data.table)

# Load common reference-based parameters
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_header.R")

# parameters specific for RCTD and panel (or disease)
UMI_min_sigma     <- snakemake@params[["UMI_min_sigma"]]
class_level       <- snakemake@params[["class_level"]] # optional, should be NULL if not provided
cores             <- snakemake@params[["cores"]]
CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]

###### end of snakemake params  ###### 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'counts' 
test_layer <- 'counts' 

# Load reference and query data
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_load_data.R") # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

# Generate reference object
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_generate_reference_obj.R")
xe_chrom_common_genes <- intersect(rownames(xe), rownames(chrom))


# Create RCTD-specific reference object 
ref_labels <-  chrom@meta.data %>% pull(annotation_level) %>% as.vector() %>% as.factor()
names(ref_labels) <- colnames(chrom)

ref.obj <- Reference(
  GetAssayData(chrom, assay = ref_assay, layer = ref_layer)[xe_chrom_common_genes, ], 
  cell_types = ref_labels, 
  min_UMI = REF_MIN_UMI, 
  require_int = !(xe@misc$sample_metadata[["segmentation_method"]] %in% c("proseg", "segger")))

# Create query object
coords   <- xe@meta.data %>% select(ST_1, ST_2)
test.obj <- SpatialRNA(
  coords,
  GetAssayData(xe, assay = xe_assay, layer = test_layer),
  require_int = snakemake@params[["is_int_counts"]]
)

# Create `class_df` if a valid annotation level is provided
class_df <- generate_class_df(
  chrom = chrom,
  annotation_level = annotation_level, 
  class_level = class_level
)


# Annotate with RCTD
# run RCTD with many cores #TODO: make it running in chuncks for large samples and gathering results afterwards
RCTD <- create.RCTD(
  test.obj,
  ref.obj,
  UMI_min = XE_MIN_UMI, #10
  counts_MIN = XE_MIN_counts, #10
  UMI_min_sigma = UMI_min_sigma, # 1, but 300 by default
  max_cores = cores,
  CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
  class_df = class_df
)

message(paste("N = ", ncol(xe) - RCTD@spatialRNA@counts@Dim[[2]], "cells were removed by create.RCTD()" ))
message(paste("N = ", nrow(xe) - RCTD@spatialRNA@counts@Dim[[1]], "genes were removed by create.RCTD()" ))

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# Get weights
weights <- as.data.frame(RCTD@results$weights)[colnames(xe), ]
weights$cell_id  <- colnames(xe)
weights          <- weights %>% select(cell_id, everything())


# Get labels 
annotations.df  <- RCTD@results$results_df[colnames(xe),]
if("scond_type" %in% colnames(annotations.df)){
  annotations.df <- annotations.df %>% mutate(second_type = scond_type)
  annotations.df$scond_type <- NULL
}
labels          <- annotations.df %>% select(first_type, second_type) 
labels$cell_id  <- colnames(xe)
labels          <- labels %>% select(cell_id, everything())

# Update second type for highly confident cells 
# Define minimum weight threshold and calculate boolean weights and candidate counts
weights_bool <- weights > 0.01
n_candidates <- rowSums(weights_bool)
high_confidence_cells <- n_candidates < 2

# Update labels based on the threshold
labels <- labels %>%
  mutate(second_type_updated = if_else(high_confidence_cells, NA_character_, second_type))

# Save annotation
saveRDS(RCTD, snakemake@output[["rds_output"]]) 
write_parquet(labels, snakemake@output[["labels"]])
write_parquet(weights, snakemake@output[["scores"]])

# Convert and Save rctd output in Py-compatible format 

# save cell_id at least in one .parquet object
results_df <- RCTD@results$results_df
results_df$cell_id <- rownames(results_df)
results_df <- results_df %>% select(cell_id, everything())
# transform score_mat into a long data.frame
score_mat <- convert_score_mat_to_long_df(RCTD@results$score_mat)
# transform singlet_scores into a long data.frame
singlet_scores <- convert_singlet_scores_to_long_df(RCTD@results$singlet_scores)

write_parquet(results_df, snakemake@output[["out_res_df"]])
write_parquet(RCTD@results$weights %>% as.data.frame(), snakemake@output[["out_w"]])
write_parquet(RCTD@results$weights_doublet %>% as.data.frame(), snakemake@output[["out_wd"]])
write_parquet(singlet_scores, snakemake@output[["out_sc"]])
write_parquet(score_mat, snakemake@output[["out_sm"]])



