# Annotate Xenium from Chromium using RCTD

library(Seurat)
library(dplyr)
library(spacexr)
source(snakemake@source("../../../scripts/_celltype_annotation/annotation_utils.R"))

###### snakemake params  ###### 

reference_path    <- snakemake@input[["reference"]] 
query_path        <- snakemake@input[["query"]] 

annotation_parameters <- get_referencebased_annotation_parameters(snakemake@params[["annotation_id"]] )
annotation_level  <- annotation_parameters[["annotation_level"]] 
reference_type    <- annotation_parameters[["reference_type"]]
annotation_mode   <- annotation_parameters[["annotation_mode"]] # So far, runs in single-cell mode only 
annotation_method <- annotation_parameters[["annotation_mode"]]

# Assays used for annotation
ref_assay         <- snakemake@params[["ref_default_assay"]] 
xe_assay          <- snakemake@params[["xe_default_assay"]] 

# Reference and query cell filtering params 
REF_MIN_UMI       <- snakemake@params[["REF_MIN_UMI"]] 
REF_MAX_UMI       <- snakemake@params[["REF_MAX_UMI"]] 
XE_MIN_UMI        <- snakemake@params[["XE_MIN_UMI"]] # mb do not provide option to specify it, but compute from Xenium object, that has been QCed already, so all cells have to be annotated?
XE_MIN_counts     <- snakemake@params[["XE_MIN_counts"]] # same as above

# up to here, all the parameters are shared across all methods and diseases

# parameters specific for RCTD and panel (or disease)
CELL_MIN_INSTANCE <- snakemake@params[["CELL_MIN_INSTANCE"]] # ideally, should de decreased for the breast panel, but lets keep it like this 
UMI_min_sigma     <- snakemake@params[["UMI_min_sigma"]] # only used in RCTD (should be passed to snakemake via _other_options?)
class_level       <- snakemake@params[["class_level"]] # optional, should be NULL if not provided

###### end of snakemake params  ###### 


# Annotation method specific, should not be modified, so not in snakemake
ref_layer  <- 'counts' 
test_layer <- 'counts' 


## Load Chromium (reference) data 
chrom <- readRDS(reference_path) # `chrom_file_path` from snakemake

## Load Xenium (query) data 
xe <- readRDS(query_path)


## Generate reference object 
chrom <- generate_reference_obj(
    chrom = chrom,
    query_features = rownames(xe),
    donor_id = xe@misc$sample_metadata[["donor"]],
    reference_type = reference_type,
    annotation_level = annotation_level,
    ref_assay = ref_assay,
    REF_MIN_UMI = REF_MIN_UMI,
    REF_MAX_UMI = REF_MAX_UMI,
    CELL_MIN_INSTANCE = CELL_MIN_INSTANCE
)

# Create reference object 
ref.obj <- Reference(GetAssayData(chrom_i, assay = ref_assay, layer = ref_layer), ref_labels, min_UMI = MIN_UMIs, require_int = (xe@misc$sample_metadata[["segmentation_method"]]!="proseg"))

# Create query object
coords   <- xe@meta.data %>% select(ST_1, ST_2)
test.obj <- SpatialRNA(coords, GetAssayData(xe, assay = xe_assay, layer = test_layer))

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
  UMI_min_sigma = UMI_min_sigma, #100, but 300 by default
  max_cores = parallel::detectCores(),
  CELL_MIN_INSTANCE = CELL_MIN_INSTANCE
)

message(paste("N = ", ncol(xe) - RCTD@spatialRNA@counts@Dim[[2]], "cells were removed by create.RCTD()" ))
message(paste("N = ", nrow(xe) - RCTD@spatialRNA@counts@Dim[[1]], "genes were removed by create.RCTD()" ))

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# Get scores, here weights 
scores <- as.data.frame(RCTD@results$weights)[colnames(xe), ]
rownames(scores) <- colnames(xe)

# Get labels 
annotations.df  <- RCTD@results$results_df[colnames(xe),]
labels          <- annotations.df %>% select(first_type, second_type) # #TODO: second_type_updated 
rownames(labels)<- colnames(xe)

# Update second type for highly confident cells 
# Define minimum weight threshold and calculate boolean scores and candidate counts
scores_bool <- scores > 0.01
n_candidates <- rowSums(scores_bool)
hc_idx <- n_candidates < 2

# Update labels based on the threshold
labels <- labels %>%
  mutate(second_type_updated = if_else(hc_idx, NA_character_, second_type))

# Save annotation
saveRDS(RCTD, skakemake@output[[1]])
write.csv(labels, skakemake@output[[2]])
write.csv(scores, skakemake@output[[3]])


