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
sample_id <- snakemake@params[["sample_id"]]
segmentation_id <- snakemake@params[["segmentation_id"]]

split_sample_id     <- strsplit(sample_id, split = "/")[[1]]
segmentation_method <- strsplit(segmentation_id, "_")[[1]][1]

xe@misc$sample_metadata <- list(
  disease = split_sample_id[1],
  gene_panel = split_sample_id[2],
  donor = split_sample_id[3],
  sample = split_sample_id[4],
  sample_id = sample_id,
  segmentation_method = segmentation_method,
  segmentation_id = segmentation_id
)

xe <- AddMetaData(xe, xe@misc$sample_metadata)
xe <- AddMetaData(xe, colnames(xe), col.name = "cell_id")
