library(arrow)
library(dplyr)

snakemake@source("../../../scripts/segmentation/_grid/grid_utils.R")

cellb_path <- file.path(snakemake@input[["xenium_bundle"]], "cell_boundaries.parquet")  # adjust if nested

message("Computing bbox from xenium bundle...\n")

bbox <- arrow::read_parquet(cellb_path, as_data_frame = FALSE) %>%
  dplyr::summarise(
    xmin = min(vertex_x, na.rm = TRUE),
    xmax = max(vertex_x, na.rm = TRUE),
    ymin = min(vertex_y, na.rm = TRUE),
    ymax = max(vertex_y, na.rm = TRUE)
  ) %>%
  collect() %>%
  unlist(use.names = TRUE)

message("Computing bbox from xenium bundle... Done\n")

message("Making Visium grid over the bbox...\n")

visium_grid <- make_visium_grid_sf(
    bbox,
    pitch_um = 100,
    spot_diameter_um = snakemake@input[["diameter"]],
    crs = NA
)
message("Making Visium grid over the bbox... Done\n")

message("Writing new cell boundaries...\n")

write_xenium_cells_geojson(visium_grid, snakemake@output[["cells"]])

message("Writing new cell boundaries... Done\n")



