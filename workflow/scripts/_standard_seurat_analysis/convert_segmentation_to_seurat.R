log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)

xr_dir <- snakemake@input[["xr_dir"]]
spatial_dimname <- snakemake@param[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = file.path(xr_dir, "outs")
)

# Create new Dimred from spatial coordinates for easier exploration, usage and ploting spatial data
st <- xe@images[["fov"]]@boundaries[["centroids"]]@coords
colnames(st) <- paste0("ST_", 1:ncol(st))
rownames(st) <- colnames(xe)

# Create Dimred from spatial coordinates
if (!is.null(coords <- xe@images[["fov"]]@boundaries[["centroids"]]@coords)) {
  if (is.matrix(coords) || is.data.frame(coords)) {
    colnames(coords) <- paste0("ST_", seq_len(ncol(coords)))
    rownames(coords) <- colnames(xe) 
    xe[[spatial_dimname]] <- CreateDimReducObject(embeddings = as.matrix(st), assay = DefaultAssay(xe))
  } else warning("Spatial coordinates are not in the expected format.")
} else warning("Spatial coordinates not found in the 'xe' object.")

  
# Save object 
saveRDS(
  xe, 
  file = file.path(snakemake@output[["xe"]])
)
