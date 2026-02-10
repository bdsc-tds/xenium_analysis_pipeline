log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat) 
library(dplyr)
library(arrow)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

default_assay <- snakemake@params[["default_assay"]]
n_dims <- snakemake@params[["n_dims"]]
resolution <- snakemake@params[["resolution"]]

dims <- 1:n_dims

# Optional sketch params (set defaults if not provided)
enable_sketching <- snakemake@params[["enable_sketching"]]
sketch_ncells <- snakemake@params[["sketch_ncells"]]
sketch_method <- snakemake@params[["sketch_method"]]
sketch_threshold <- snakemake@params[["sketch_threshold"]]

# Read post QC Seurat
xe <- readRDS(snakemake@input[[1]])

ncells <- ncol(xe)

if (enable_sketching & ncells > sketch_threshold) {
  message("Large object detected (", ncells, " cells). Running Seurat sketch workflow.")
  
  # 1) Build sketch assay (stores a representative subset of cells)
  xe <- SketchData(
    object = xe,
    ncells = sketch_ncells,
    method = sketch_method,
    assay = default_assay,
    new.assay.name = "sketch"
  )
  
  # Reference = sketch cells only
  sketch_cells <- Cells(xe, assay = "sketch")
  #ref <- subset(xe, cells = sketch_cells)
  
  DefaultAssay(xe) <- "sketch"
  
  # 2) Run full graph-based workflow on the sketch only
  xe <- xe |>
    ScaleData() |>
    RunPCA(
      npcs = n_dims,
      reduction.name = "sketching_pca",
      reduction.key = "skPC_"
    ) |>
    FindNeighbors(
      reduction = "sketching_pca",
      dims = dims
    ) |>
    FindClusters(
      resolution = resolution
    ) |>
    RunUMAP(
      reduction = "sketching_pca",
      dims = dims,
      return.model = TRUE,
      reduction.name = "sketching_umap",
      reduction.key = "skUMAP_"
    )
  
  # 3) Project labels + UMAP back to all cells (query = full object)
  DefaultAssay(xe) <- default_assay
  
  xe <- ProjectData(
    object = xe,
    assay = default_assay,
    full.reduction = "pca",
    sketched.assay = "sketch",
    sketched.reduction = "sketching_pca",
    umap.model = "sketching_umap",
    dims = dims,
    refdata = list(clusters = "seurat_clusters")
  )
  
  xe@reductions[["umap"]] <- xe@reductions[["full.sketching_umap"]]
  colnames(xe@reductions$umap@cell.embeddings) <- c("umap_1", "umap_2")
  xe@reductions[["umap"]]@key <- "umap_"
  
  xe@reductions[["pca"]]@stdev <- xe@reductions[["sketching_pca"]]@stdev
  xe@reductions[["pca"]]@feature.loadings <- xe@reductions[["sketching_pca"]]@feature.loadings
  
  # clean up
  xe[["sketch"]] <- NULL
  xe@reductions[["sketching_pca"]] <- NULL
  xe@reductions[["sketching_umap"]] <- NULL
  xe@reductions[["full.sketching_umap"]] <- NULL
  
  # Save some sketch-specific metadata too (useful for debugging)
  xe@misc$sketch_meta <- list(
    sketch_ncells = sketch_ncells,
    sketch_method = sketch_method,
    sketch_cells = sketch_cells
  )
  
} else {
  message("Object size (", ncells, " cells) <= ", sketch_threshold, " or sketching is disabled. Running standard workflow.")
  
  xe <- xe |>
    RunPCA(npcs = n_dims) |>
    FindNeighbors(dims = dims) |>
    FindClusters(resolution = resolution) |>
    RunUMAP(dims = dims)
  
  xe@meta.data <- xe@meta.data |>
    mutate(clusters = seurat_clusters)
}

# Save parameters in seurat object
xe@misc$standard_seurat_analysis_meta <- c(
  xe@misc$standard_seurat_analysis_meta,
  list(
    n_dims = n_dims,
    dims = dims,
    resolution = resolution
  )
)

# Save results
saveRDS(
  xe, 
  file = file.path(snakemake@output[["obj"]])
)

write_parquet(
  data.frame(
    cell = colnames(xe)
  ),
  sink = file.path(snakemake@output[["cells"]])
)

write_parquet(
  as.data.frame(
    Embeddings(xe, reduction = "pca")
  ),
  sink = file.path(snakemake@output[["pca"]])
)

write_parquet(
  as.data.frame(
    Embeddings(xe, reduction = "umap")
  ),
  sink = file.path(snakemake@output[["umap"]])
)
