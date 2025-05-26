# Add sample-related metadata to Seurat object
xe@misc$sample_metadata <- list(
  segmentation_id = snakemake@params[["segmentation_id"]],
  segmentation_method = snakemake@params[["segmentation_method"]],
  sample_id = snakemake@params[["sample_id"]],
  condition = snakemake@params[["condition"]],
  gene_panel = snakemake@params[["gene_panel"]],
  donor = snakemake@params[["donor"]],
  sample = snakemake@params[["sample"]],
  normalisation_id = snakemake@params[["normalisation_id"]],
  annotation_id = snakemake@params[["annotation_id"]],
  annotation_approach = snakemake@params[["annotation_approach"]],
  annotation_reference_name = snakemake@params[["annotation_reference_name"]],
  annotation_method = snakemake@params[["annotation_method"]],
  annotation_level = snakemake@params[["annotation_level"]],
  annotation_mode = snakemake@params[["annotation_mode"]],
  count_correction_id = snakemake@params[["count_correction_id"]]
)

xe <- AddMetaData(xe, xe@misc$sample_metadata)
