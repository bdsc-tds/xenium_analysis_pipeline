# Generating reference data for the reference-based cell-type annotation

# Validate metadata
if (!"donor" %in% names(xe@misc$sample_metadata)) stop("The 'donor' column is missing in Xenium metadata!")
if (!annotation_level %in% colnames(chrom@meta.data)) stop(paste("Annotation level", annotation_level, "is missing in reference metadata!"))

message("Generating reference object...")
chrom <- generate_reference_obj(
  chrom = chrom,
  query_features = rownames(xe),
  donor_id = xe@misc$sample_metadata[["donor"]],
  # reference_type = reference_type,
  annotation_level = annotation_level,
  ref_assay = ref_assay,
  MIN_UMI_ref = REF_MIN_UMI,
  MAX_UMI_ref = REF_MAX_UMI,
  CELL_MIN_INSTANCE = CELL_MIN_INSTANCE
)
message("Reference object generated...")