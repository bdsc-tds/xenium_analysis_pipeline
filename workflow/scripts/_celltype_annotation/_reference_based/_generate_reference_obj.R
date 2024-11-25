# Generating reference data for the reference-based cell-type annotation

# Validate metadata
if (!"donor" %in% colnames(xe@misc$sample_metadata)) stop("The 'donor' column is missing in Xenium metadata!")
if (!annotation_level %in% colnames(chrom@meta.data)) stop(paste("Annotation level", annotation_level, "is missing in reference metadata!"))

message("Generating reference object...")
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