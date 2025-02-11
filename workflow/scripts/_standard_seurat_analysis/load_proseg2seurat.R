log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

snakemake@source("../../scripts/utils/run_time_utils.R")

library(arrow)
library(SeuratObject)
library(Seurat)
library(Matrix)
library(dplyr)

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = snakemake@params[["data_dir"]]
)

mapping <- read_parquet(snakemake@input[["mapping"]])

# Adjust indexing from Python 0-based to R 1-based
mapping[[snakemake@params[["proseg_cell_id_col_name"]]]] <- mapping[[snakemake@params[["proseg_cell_id_col_name"]]]] + 1L

if (snakemake@params[["use_mode_counts"]]) {
  if (snakemake@params[["use_mapping"]]) {
    xe <- xe[, mapping[[snakemake@params[["xr_cell_id_col_name"]]]]]
  }
} else {
  # Adjust indexing from Python 0-based to R 1-based
  cell_metadata <- read.csv(snakemake@input[["cell_metadata"]]) %>%
    mutate(cell = cell + 1) %>%
    arrange(cell)

  expected_counts <- read.csv(snakemake@input[["expected_counts"]])

  if (snakemake@params[["use_mapping"]]) {
    cell_metadata <- cell_metadata[cell_metadata$cell %in% mapping[[snakemake@params[["proseg_cell_id_col_name"]]]], ]
    expected_counts <- expected_counts[mapping[[snakemake@params[["proseg_cell_id_col_name"]]]], ]
  }

  if (nrow(cell_metadata) != nrow(expected_counts)) {
    stop("Error! Number of cells in cell metadata and expected counts in the results of Proseg do not match.")
  }

  set.seed(42)
  new_cell_ids <- generate_xr_style_cell_names(nrow(expected_counts))
  gene_ids <- colnames(expected_counts)

  expected_counts <- t(as.matrix(expected_counts))

  # Create a new Seurat object for expected counts
  new_xe <- CreateSeuratObject(
    counts = expected_counts,
    assay = "Xenium"
  )
  colnames(new_xe) <- new_cell_ids
  rownames(new_xe) <- gene_ids

  coords <- CreateCentroids(
    coords = cell_metadata[, c("centroid_x", "centroid_y")],
    nsides = xe@images$fov@boundaries$centroids@nsides,
    radius = xe@images$fov@boundaries$centroids@radius,
    theta = xe@images$fov@boundaries$centroids@theta
  )
  coords <- RenameCells(coords, new_cell_ids)

  new_xe@images <- list(
    fov = CreateFOV(
      coords,
      molecules = xe@images$fov@molecules[[1]],
      assay = xe@images$fov@assay,
      key = xe@images$fov@key
    )
  )

  missing_assays <- Assays(xe)[!Assays(xe) %in% Assays(new_xe)]
  names(missing_assays) <- missing_assays

  new_assays <- sapply(
    missing_assays,
    function(x) {
      assay <- GetAssayData(xe, assay = x)

      CreateAssayObject(
        counts = Matrix(
            0,
            ncol = ncol(new_xe),
            nrow = nrow(assay),
            dimnames = list(
                rownames(assay),
                new_cell_ids
            )
        )
      )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  for (i in names(new_assays)) {
    new_xe[[i]] <- new_assays[[i]]
  }

  xe <- new_xe
}

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
