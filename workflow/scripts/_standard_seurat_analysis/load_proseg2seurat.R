log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

snakemake@source("../../scripts/utils/run_time_utils.R")

library(arrow)

spatial_dimname <- snakemake@params[["spatial_dimname"]]

xe <- LoadXenium(
  data.dir = snakemake@params[["input_data"]][["dir"]],
  molecule.coordinates = FALSE
)

if (grepl(".+\\.csv.*", snakemake@params[["input_data"]][["mapping"]])) {
  mapping <- read.csv(snakemake@params[["input_data"]][["mapping"]], check.names=FALSE)
  mapping[["proseg_cell_id_bak"]] <- mapping[[snakemake@params[["proseg_cell_id_col_name"]]]]
  mapping[[snakemake@params[["proseg_cell_id_col_name"]]]] <- as.integer(gsub(
    "^cell-",
    "",
    mapping[[snakemake@params[["proseg_cell_id_col_name"]]]]
  ))
} else if (grepl(".+\\.parquet$", snakemake@params[["input_data"]][["mapping"]])) {
  mapping <- read_parquet(snakemake@params[["input_data"]][["mapping"]])
} else {
  stop("Unsupported mapping file format. Supported formats are: '.parquet', '.csv.gz' or 'csv'.")
}

# Adjust indexing from Python 0-based to R 1-based
mapping[[snakemake@params[["proseg_cell_id_col_name"]]]] <- mapping[[snakemake@params[["proseg_cell_id_col_name"]]]] + 1L

if (snakemake@params[["use_mode_counts"]]) {
  if (snakemake@params[["use_mapping"]]) {
    xe <- xe[, mapping[[snakemake@params[["xr_cell_id_col_name"]]]]]
  }
} else {
  # Adjust indexing from Python 0-based to R 1-based
  cell_metadata <- read.csv(snakemake@params[["input_data"]][["cell_metadata"]]) %>%
    mutate(cell = cell + 1) %>%
    arrange(cell)

  expected_counts <- read.csv(snakemake@params[["input_data"]][["expected_counts"]])
  expected_counts <- expected_counts[, colnames(expected_counts)[!grepl(
    snakemake@params[["control_gene_pat"]],
    colnames(expected_counts)
  )]]

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
    cell_metadata[, c("centroid_x", "centroid_y")],
    "proseg-"
  )
}

snakemake@source("../../scripts/_standard_seurat_analysis/_post_seurat_load_xenium.R")
snakemake@source("../../scripts/_standard_seurat_analysis/_add_metadata_post_seurat_load_xenium.R")

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
