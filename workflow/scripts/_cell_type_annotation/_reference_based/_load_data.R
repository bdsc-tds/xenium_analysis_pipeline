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
### REFERENCE
# Ensure ref_assay is a character string
if (!is.character(ref_assay) || length(ref_assay) == 0) {
  stop("Error: ref_assay must be a non-empty character string.")
}

# Get assays from chrom and ensure it is a character vector
assays_list <- Assays(chrom)

if (is.null(assays_list) || length(assays_list) == 0) {
  stop("Error: No assays found in the reference object.")
}

# Check if ref_assay is present in assays_list
if (!ref_assay %in% as.character(assays_list)) {
  stop("Error: Assay '", ref_assay, "' is not found in the reference object. Please update the object and re-run.")
} else {
  DefaultAssay(chrom) <- ref_assay
}

### QUERY
# Ensure xe_assay is a non-empty character string
if (!is.character(xe_assay) || length(xe_assay) == 0) {
  stop("Error: xe_assay must be a non-empty character string.")
}

# Get assays from xe and ensure it is a character vector
assays_list <- Assays(xe)

if (is.null(assays_list) || length(assays_list) == 0) {
  stop("Error: No assays found in the xenium object.")
}

# Check if xe_assay is present in assays_list
if (!xe_assay %in% as.character(assays_list)) {
  stop("Error: Assay '", xe_assay, "' is not found in the xenium object. Please update the object and re-run.")
} else {
  DefaultAssay(xe) <- xe_assay
}

message("All required assays exist...")