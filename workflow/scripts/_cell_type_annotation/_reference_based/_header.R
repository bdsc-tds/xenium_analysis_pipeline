# header.R
# Shared parameters for reference-based cell type annotation
# Loads Snakemake inputs and parameters for downstream scripts
snakemake@source("../../../scripts/_cell_type_annotation/annotation_utils.R")

###### snakemake params  ###### 

reference_path    <- snakemake@input[["reference"]] 
query_path        <- snakemake@input[["query"]] 

annotation_parameters <- get_referencebased_annotation_parameters(snakemake@params[["annotation_id"]] )
annotation_level  <- annotation_parameters[["annotation_level"]] 
reference_type    <- annotation_parameters[["reference_type"]]
annotation_mode   <- annotation_parameters[["annotation_mode"]] # So far, runs in single-cell mode only 
annotation_method <- annotation_parameters[["annotation_method"]]

# Assays used for annotation
ref_assay         <- snakemake@params[["ref_assay"]] 
xe_assay          <- snakemake@params[["xe_assay"]] 

# Reference and query cell filtering params 
REF_MIN_UMI       <- snakemake@params[["REF_MIN_UMI"]] 
REF_MAX_UMI       <- snakemake@params[["REF_MAX_UMI"]] 
XE_MIN_UMI        <- snakemake@params[["XE_MIN_UMI"]] # mb do not provide option to specify it, but compute from Xenium object, that has been QCed already, so all cells have to be annotated?
XE_MIN_counts     <- snakemake@params[["XE_MIN_counts"]] # same as above

CELL_MIN_INSTANCE <- snakemake@params[["CELL_MIN_INSTANCE"]] # ideally, should de decreased for the breast panel, but lets keep it like this 

# up to here, all the parameters are shared across all methods and diseases
