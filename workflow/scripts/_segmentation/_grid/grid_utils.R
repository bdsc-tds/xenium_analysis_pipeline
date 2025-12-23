library(sf)

make_visium_grid_sf <- function(
    bbox,
    pitch_um = 100,
    spot_diameter_um = 55,
    crs = NA
){
  stopifnot(is.numeric(bbox), length(bbox) == 4)
  names(bbox) <- names(bbox) #%||% c("xmin","ymin","xmax","ymax")
  xmin <- bbox[["xmin"]]; ymin <- bbox[["ymin"]]; xmax <- bbox[["xmax"]]; ymax <- bbox[["ymax"]]
  if (xmin > xmax) { tmp <- xmin; xmin <- xmax; xmax <- tmp }
  if (ymin > ymax) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
  
  # Visium (classic) hex packing:
  dx <- pitch_um
  dy <- pitch_um * sqrt(3) / 2   # ~86.6025 um
  xoff <- dx / 2                 # shift every other row by 50 um
  r <- spot_diameter_um / 2
  
  # expand bbox slightly so buffers still cover the bbox edges
  xmin2 <- xmin - dx
  xmax2 <- xmax + dx
  ymin2 <- ymin - dy
  ymax2 <- ymax + dy
  
  ys <- seq(ymin2, ymax2, by = dy)
  rows <- seq_along(ys)
  
  centers <- do.call(rbind, lapply(rows, function(i) {
    y <- ys[i]
    x0 <- if (i %% 2 == 0) xmin2 + xoff else xmin2
    xs <- seq(x0, xmax2, by = dx)
    cbind(x = xs, y = rep(y, length(xs)))
  }))
  
  # keep only centers whose (circular) spot intersects the bbox
  # quick filter: center within bbox expanded by radius
  keep <- centers[, "x"] >= (xmin - r) & centers[, "x"] <= (xmax + r) &
    centers[, "y"] >= (ymin - r) & centers[, "y"] <= (ymax + r)
  centers <- centers[keep, , drop = FALSE]
  
  # sf points
  pts <- sf::st_as_sf(
    data.frame(spot_id = seq_len(nrow(centers)), centers),
    coords = c("x","y"),
    crs = crs
  )
  
  # sf spot polygons as circles (buffer in same units: µm)
  spots <- sf::st_buffer(pts, dist = r)
  
  list(centers = pts, spots = spots)
}

map_cells_to_spots <- function(
    visium_grid_sf,
    cells_df
){
  
  cells_sf <- st_as_sf(cells_df, coords = c("x", "y"), crs = st_crs(visium_grid_sf$spots))
  
  cell_to_spot <- st_join(
    cells_sf,
    visium_grid_sf$spots %>% select(spot_id),
    join = st_within,
    left = TRUE
  )
}


write_xenium_cells_geojson <- function(g, out_geojson,
                                       id_col = "spot_id",
                                       make_ids_char = TRUE,
                                       ensure_valid = TRUE) {
  stopifnot(is.list(g), "spots" %in% names(g))
  stopifnot(inherits(g$spots, "sf"))
  
  suppressPackageStartupMessages({
    library(sf)
    library(jsonlite)
  })
  
  cells <- g$spots
  
  # Ensure we have an ID column
  if (!id_col %in% names(cells)) {
    cells[[id_col]] <- seq_len(nrow(cells))
  }
  if (make_ids_char) {
    cells[[id_col]] <- as.character(cells[[id_col]])
  }
  
  # Ensure polygon geometry (buffer makes polygons already)
  if (ensure_valid) {
    cells <- sf::st_make_valid(cells)
  }
  # Cast to POLYGON (not MULTIPOLYGON) to match your example
  cells <- sf::st_cast(cells, "POLYGON", warn = FALSE)
  
  # Extract coordinates as plain R lists (no CRS info)
  geom_list <- lapply(sf::st_geometry(cells), function(gm) {
    # Polygon: list(ring_matrix, ring_matrix, ...)
    rings <- lapply(gm, function(ring) {
      # ring is matrix [n,2]
      unname(split(ring, row(ring)))  # list of c(x,y)
    })
    list(type = "Polygon", coordinates = rings)
  })
  
  features <- lapply(seq_len(nrow(cells)), function(i) {
    list(
      type = "Feature",
      id = cells[[id_col]][i],
      geometry = geom_list[[i]],
      # If you DO want properties, uncomment:
      properties = list(objectType = "cell")
    )
  })
  
  fc <- list(type = "FeatureCollection", features = features)
  jsonlite::write_json(fc, out_geojson, auto_unbox = TRUE, pretty = FALSE)
  invisible(out_geojson)
}

