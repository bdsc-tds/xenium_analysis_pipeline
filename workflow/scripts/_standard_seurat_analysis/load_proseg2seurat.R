log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

snakemake@source("../../scripts/utils/run_time_utils.R")

library(arrow)

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = snakemake@params[["data_dir"]],
  molecule.coordinates = FALSE
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

  xe <- replace_counts_in_seurat(
    xe,
    expected_counts,
    cell_metadata[, c("centroid_x", "centroid_y")]
  )
}

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
