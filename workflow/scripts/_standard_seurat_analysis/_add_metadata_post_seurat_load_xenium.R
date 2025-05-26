# Add sample-related metadata to Seurat object
xe@misc$sample_metadata <- list(
  segmentation_id = snakemake@params[["segmentation_id"]],
  segmentation_method = snakemake@params[["segmentation_method"]],
  sample_id = snakemake@params[["sample_id"]],
  condition = snakemake@params[["condition"]],
  gene_panel = snakemake@params[["gene_panel"]],
  donor = snakemake@params[["donor"]],
  sample = snakemake@params[["sample"]]
)

xe <- AddMetaData(xe, xe@misc$sample_metadata)
