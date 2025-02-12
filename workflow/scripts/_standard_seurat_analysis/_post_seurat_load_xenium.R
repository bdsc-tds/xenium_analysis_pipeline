# Create Dimred from spatial coordinates
if (!is.null(coords <- xe@images[["fov"]]@boundaries[["centroids"]]@coords)) {
  if (is.matrix(coords) || is.data.frame(coords)) {
    colnames(coords) <- paste0("ST_", seq_len(ncol(coords)))
    rownames(coords) <- colnames(xe) 
    xe[[spatial_dimname]] <- CreateDimReducObject(embeddings = as.matrix(coords), assay = DefaultAssay(xe))
    
    # Add spatial coordinates into metadata 
    xe <- AddMetaData(xe, coords)
  } else warning("Spatial coordinates are not in the expected format.")
} else warning("Spatial coordinates not found in the Seurat object.")

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
xe <- AddMetaData(xe, colnames(xe), col.name = "cell_id")
