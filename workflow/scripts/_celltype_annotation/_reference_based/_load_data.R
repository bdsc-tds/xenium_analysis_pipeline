# Loading of the reference and query data for the reference-based cell-type annotation

# Validate input files
if (!file.exists(reference_path)) stop(paste("Reference file not found:", reference_path))
if (!file.exists(query_path)) stop(paste("Query file not found:", query_path))

message("Loading data...")
chrom <- readRDS(reference_path)
xe <- readRDS(query_path)
message("Data loaded.")

