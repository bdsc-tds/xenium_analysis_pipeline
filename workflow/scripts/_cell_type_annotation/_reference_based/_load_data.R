# Loading of the reference and query data for the reference-based cell-type annotation

# Validate input files
if (!file.exists(reference_path)) stop(paste("Reference file not found:", reference_path))
if (!file.exists(query_path)) stop(paste("Query file not found:", query_path))

message("Loading data...")
chrom <- readRDS(reference_path)
xe <- readRDS(query_path)
message("Data loaded.")

message("Checking required assays...")
# Make sure required assays exist
if(!ref_assay %in% Assays(chrom)){
  stop("Assay", ref_assay, "is not found in the reference object, please update the object and re-run!")
} else {
  DefaultAssay(chrom) <- ref_assay
}

if(!xe_assay %in% Assays(xe)){
  stop("Assay", xe_assay, "is not found in the xenium object, please update the object and re-run!")
} else {
  DefaultAssay(xe)    <- xe_assay
}
message("All required assays exist...")