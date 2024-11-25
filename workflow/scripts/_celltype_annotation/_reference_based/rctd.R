# Annotate Xenium from Chromium using RCTD

library(Seurat)
library(dplyr)
library(spacexr)

# Load common reference-based parameters
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/header.R")) 

# parameters specific for RCTD and panel (or disease)
CELL_MIN_INSTANCE <- snakemake@params[["CELL_MIN_INSTANCE"]] # ideally, should de decreased for the breast panel, but lets keep it like this 
UMI_min_sigma     <- snakemake@params[["UMI_min_sigma"]] # only used in RCTD (should be passed to snakemake via _other_options?)
class_level       <- snakemake@params[["class_level"]] # optional, should be NULL if not provided

###### end of snakemake params  ###### 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'counts' 
test_layer <- 'counts' 

# Load reference and query data
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_load_data.R")) # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

# Generate reference object
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_generate_reference_obj.R")) 

# Create RCTD-specific reference object 
ref.obj <- Reference(
  GetAssayData(chrom, assay = ref_assay, layer = ref_layer), 
  cell_types = chrom@meta.data %>% pull(annotation_level) %>% as.vector() %>% as.factor(), 
  min_UMI = REF_MIN_UMI, 
  require_int = (xe@misc$sample_metadata[["segmentation_method"]]!="proseg"))

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
high_confidence_cells <- n_candidates < 2

# Update labels based on the threshold
labels <- labels %>%
  mutate(second_type_updated = if_else(high_confidence_cells, NA_character_, second_type))

# Save annotation
saveRDS(RCTD, snakemake@output[[1]])
write.csv(labels, snakemake@output[[2]])
write.csv(scores, snakemake@output[[3]])


