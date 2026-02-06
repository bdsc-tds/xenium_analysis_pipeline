log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(arrow)
library(dplyr)

snakemake@source(file.path(snakemake@scriptdir, "grid_utils.R"))

cellb_path <- file.path(snakemake@input[["xenium_bundle"]], "cell_boundaries.parquet")  # adjust if nested

message("Xenium bundle path: ", cellb_path)

if(!file.exists(cellb_path)){
  stop("Xenium bundle ", cellb_path, "does not exist! \n")
}

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

message("Making VisiumHD grid over the bbox...\n")

visiumHD_grid <- make_visiumhd_grid_sf(
    bbox,
    bin_size_um = snakemake@params[["bin_size"]],
    crs = NA
)
message("Making VisiumHD grid over the bbox... Done\n")

message("Writing new cell boundaries...\n")

write_xenium_cells_geojson(visiumHD_grid, snakemake@output[["cells"]])

message("Writing new cell boundaries... Done\n")



