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


#' Create a Visium HD square-bin grid as `sf` objects
#'
#' Generates a square (axis-aligned) grid of Visium HD-style bins covering a
#' bounding box in micron units. Returns both bin centers (POINT) and bin
#' geometries (POLYGON), suitable for spatial joins (e.g., mapping Xenium cells
#' to bins with `sf::st_join()`).
#'
#' The grid is defined by a bin size (`bin_size_um`) and an anchor `origin` that
#' controls how bin edges align to the bounding box. Optionally, the bounding
#' box can be expanded by one bin to ensure edge bins intersect the bbox.
#'
#' @param bbox Numeric vector of length 4 with names `xmin`, `ymin`, `xmax`,
#'   `ymax` defining the bounding box (in microns).
#' @param bin_size_um Positive numeric scalar. Side length of each square bin
#'   in microns.
#' @param crs Coordinate reference system passed to `sf::st_as_sf()`. Use `NA`
#'   for planar (unitless) coordinates; for real-world CRS supply an EPSG code
#'   or `sf::st_crs()` object.
#' @param origin Character scalar controlling the anchor of bin edges. Either
#'   `"xmin"` (default) or `"xmax"` for the x-axis, and `"ymin"` (default) or
#'   `"ymax"` for the y-axis. Values are combined as `c("xmin", "ymin")`,
#'   `c("xmax", "ymin")`, etc. Internally only `"xmin"`/`"xmax"` and
#'   `"ymin"`/`"ymax"` are used; the default aligns bins from the lower-left.
#' @param expand Logical. If `TRUE` (default), expands the bbox by one bin in
#'   all directions before generating the grid to improve edge coverage.
#'
#' @return A named list with two `sf` objects:
#' \describe{
#'   \item{centers}{`sf` POINT features with columns `spot_id` and center
#'     coordinates.}
#'   \item{spots}{`sf` POLYGON features (square bins) with column `spot_id`.}
#' }
#'
#' @details
#' Filtering is performed to keep only bins whose squares intersect the original
#' bbox. This is approximated by retaining centers within the bbox expanded by
#' half a bin in each direction.
#'
#' @examples
#' bbox <- c(xmin = 0, ymin = 0, xmax = 5000, ymax = 4000)
#'
#' g8 <- make_visiumhd_grid_sf(bbox, bin_size_um = 8)
#' g16 <- make_visiumhd_grid_sf(bbox, bin_size_um = 16)
#'
#' # Plot polygons
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   plot(sf::st_geometry(g8$spots), border = "grey70")
#'   plot(sf::st_geometry(g8$centers), add = TRUE, pch = 16, cex = 0.3)
#' }
#'
make_visiumhd_grid_sf <- function(
    bbox,
    bin_size_um = 8,
    crs = NA,
    origin = c("xmin", "ymin"),
    expand = TRUE
) {
  stopifnot(is.numeric(bbox), length(bbox) == 4)
  names(bbox) <- names(bbox)
  xmin <- bbox[["xmin"]]; ymin <- bbox[["ymin"]]
  xmax <- bbox[["xmax"]]; ymax <- bbox[["ymax"]]
  
  if (xmin > xmax) { tmp <- xmin; xmin <- xmax; xmax <- tmp }
  if (ymin > ymax) { tmp <- ymin; ymin <- ymax; ymax <- tmp }
  
  stopifnot(is.numeric(bin_size_um), length(bin_size_um) == 1, bin_size_um > 0)
  
  # Optionally expand bbox by one bin so edge bins still cover bbox boundary
  if (expand) {
    xmin2 <- xmin - bin_size_um
    xmax2 <- xmax + bin_size_um
    ymin2 <- ymin - bin_size_um
    ymax2 <- ymax + bin_size_um
  } else {
    xmin2 <- xmin
    xmax2 <- xmax
    ymin2 <- ymin
    ymax2 <- ymax
  }
  
  # Define anchor/origin for the grid
  origin <- match.arg(origin)
  x_anchor <- if (origin == "xmin") xmin2 else xmax2
  y_anchor <- if (origin == "ymin") ymin2 else ymax2
  
  # Create sequences of BIN EDGES, then derive centers
  # We build the grid so that bin edges are aligned to the anchor.
  if (origin == "xmin") {
    x_edges <- seq(x_anchor, xmax2 + bin_size_um, by = bin_size_um)
  } else {
    x_edges <- seq(x_anchor, xmin2 - bin_size_um, by = -bin_size_um)
    x_edges <- sort(x_edges)
  }
  
  if (origin == "ymin") {
    y_edges <- seq(y_anchor, ymax2 + bin_size_um, by = bin_size_um)
  } else {
    y_edges <- seq(y_anchor, ymin2 - bin_size_um, by = -bin_size_um)
    y_edges <- sort(y_edges)
  }
  
  # Centers are midpoints between edges
  x_centers <- (x_edges[-length(x_edges)] + x_edges[-1]) / 2
  y_centers <- (y_edges[-length(y_edges)] + y_edges[-1]) / 2
  
  centers <- as.matrix(expand.grid(x = x_centers, y = y_centers))
  
  # Keep only bins whose SQUARE intersects bbox:
  # Quick filter using center within bbox expanded by half-bin
  half <- bin_size_um / 2
  keep <- centers[, "x"] >= (xmin - half) & centers[, "x"] <= (xmax + half) &
    centers[, "y"] >= (ymin - half) & centers[, "y"] <= (ymax + half)
  centers <- centers[keep, , drop = FALSE]
  
  # sf points (centers)
  pts <- sf::st_as_sf(
    data.frame(spot_id = seq_len(nrow(centers)), centers),
    coords = c("x", "y"),
    crs = crs
  )
  
  # Build square polygons from centers (axis-aligned)
  # Each square: (x±half, y±half)
  sq <- lapply(seq_len(nrow(centers)), function(i) {
    x <- centers[i, "x"]
    y <- centers[i, "y"]
    ring <- matrix(
      c(
        x - half, y - half,
        x + half, y - half,
        x + half, y + half,
        x - half, y + half,
        x - half, y - half
      ),
      ncol = 2,
      byrow = TRUE
    )
    sf::st_polygon(list(ring))
  })
  
  spots <- sf::st_sf(
    spot_id = pts$spot_id,
    geometry = sf::st_sfc(sq, crs = sf::st_crs(pts))
  )
  
  list(centers = pts, spots = spots)
}
