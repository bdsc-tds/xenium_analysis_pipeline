log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Turn the rctd-py results into the spacexr RCTD object and the files rctd.R writes, so
# runPostprocessRCTD and every downstream rule read them unchanged.
# The object is built with create.RCTD (same query, reference cells and settings as the Python run),
# then its results are filled with the rctd-py fit.

library(Matrix)
library(dplyr)
library(spacexr)
library(arrow)
library(rhdf5)
library(jsonlite)

# Read a 10x v3 HDF5 count matrix (genes x cells)
read_10x_h5 <- function(path) {
  m <- h5read(path, "matrix", bit64conversion = "double")
  counts <- sparseMatrix(
    i    = as.integer(m$indices) + 1L,
    p    = as.integer(m$indptr),
    x    = as.numeric(m$data),
    dims = as.integer(m$shape)
  )
  rownames(counts) <- as.character(m$features$name)
  colnames(counts) <- as.character(m$barcodes)
  counts
}

annotation_level  <- snakemake@params[["annotation_level"]]
UMI_min_sigma     <- snakemake@params[["UMI_min_sigma"]]
CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]
REF_MIN_UMI       <- snakemake@params[["REF_MIN_UMI"]]
XE_MIN_UMI        <- snakemake@params[["XE_MIN_UMI"]]
XE_MIN_counts     <- snakemake@params[["XE_MIN_counts"]]
is_int_counts     <- snakemake@params[["is_int_counts"]]

run_info <- fromJSON(snakemake@input[["run_info"]])
message(sprintf("rctd-py %s (torch %s, %s), sigma %d", run_info$rctd_py_version, run_info$torch_version, run_info$device, run_info$sigma))

# Query: counts and coordinates, as SpatialRNA in rctd.R
xe_counts <- read_10x_h5(snakemake@input[["query_counts"]])
xe_meta   <- as.data.frame(read_parquet(snakemake@input[["query_metadata"]]))
rownames(xe_meta) <- xe_meta$cell_id
xe_cells  <- colnames(xe_counts)
test.obj  <- SpatialRNA(
  xe_meta[xe_cells, c("ST_1", "ST_2")],
  xe_counts,
  require_int = is_int_counts
)

# Reference: the cells and labels the Python run kept, on the genes shared with the query
ref_cells   <- as.data.frame(read_parquet(snakemake@input[["reference_cells"]]))
ref_counts  <- read_10x_h5(snakemake@input[["ref_counts"]])
common_genes <- intersect(rownames(xe_counts), rownames(ref_counts))
ref_labels  <- as.factor(ref_cells$cell_type)
names(ref_labels) <- ref_cells$cell_id
ref.obj <- Reference(
  ref_counts[common_genes, ref_cells$cell_id],
  cell_types  = ref_labels,
  min_UMI     = REF_MIN_UMI,
  require_int = is_int_counts
)
rm(ref_counts)

# class_df as generate_class_df builds it (NULL without a class level)
class_tbl <- as.data.frame(read_parquet(snakemake@input[["class_df"]]))
class_df  <- NULL
if (nrow(class_tbl) > 0) {
  class_df <- data.frame(class = class_tbl$class, row.names = class_tbl$cell_type)
}

# RCTD object with the settings of rctd.R (cell type info, gene lists, filtered query)
RCTD <- create.RCTD(
  test.obj,
  ref.obj,
  UMI_min           = XE_MIN_UMI,
  counts_MIN        = XE_MIN_counts,
  UMI_min_sigma     = UMI_min_sigma,
  max_cores         = 1,
  CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
  class_df          = class_df
)
cell_type_names <- RCTD@cell_type_info$info[[2]]
if (!identical(sort(cell_type_names), sort(run_info$cell_type_names))) {
  stop("Cell types of the R reference and of the rctd-py run differ: ",
       paste(cell_type_names, collapse = ","), " vs ", paste(run_info$cell_type_names, collapse = ","))
}

# rctd-py results, restricted to the cells create.RCTD kept
py_res <- as.data.frame(read_parquet(snakemake@input[["results_df"]]))
py_w   <- as.data.frame(read_parquet(snakemake@input[["weights"]]))
py_wd  <- as.data.frame(read_parquet(snakemake@input[["weights_doublet"]]))
r_cells <- colnames(RCTD@spatialRNA@counts)
fitted  <- py_res$cell_id[py_res$cell_id %in% r_cells]
message(sprintf("rctd-py fitted %d cells, create.RCTD kept %d, %d in both", nrow(py_res), length(r_cells), length(fitted)))
if (length(fitted) < 0.99 * length(r_cells)) stop("The rctd-py and spacexr cell filters disagree on more than 1% of cells.")
rownames(py_res) <- py_res$cell_id
rownames(py_w)   <- py_w$cell_id
rownames(py_wd)  <- py_wd$cell_id

# results_df with spacexr's column names and types (scond_type, character types, factor spot_class)
results_df <- data.frame(
  spot_class    = factor(py_res[fitted, "spot_class"], levels = c("reject", "singlet", "doublet_certain", "doublet_uncertain")),
  first_type    = as.character(py_res[fitted, "first_type"]),
  scond_type    = as.character(py_res[fitted, "scond_type"]),
  first_class   = py_res[fitted, "first_class"],
  second_class  = py_res[fitted, "second_class"],
  min_score     = py_res[fitted, "min_score"],
  singlet_score = py_res[fitted, "singlet_score"],
  conv_all      = NA,
  conv_doublet  = NA,
  row.names     = fitted
)

# Weights as in spacexr: cells x cell types, and the two doublet weights
weights <- Matrix(as.matrix(py_w[fitted, cell_type_names]), sparse = TRUE)
dimnames(weights) <- list(fitted, cell_type_names)
weights_doublet <- Matrix(as.matrix(py_wd[fitted, c("first_type", "second_type")]), sparse = TRUE)
dimnames(weights_doublet) <- list(fitted, c("first_type", "second_type"))

# Per-type singlet scores: one named vector per fitted cell, in results_df order (results$singlet_scores)
py_sc <- as.data.frame(read_parquet(snakemake@input[["singlet_scores"]]))
py_sc$cell_id <- py_res$cell_id[py_sc$cell_index]
py_sc <- py_sc[py_sc$cell_id %in% fitted, ]
singlet_scores <- split(setNames(py_sc$score, py_sc$cell_type), factor(py_sc$cell_id, levels = fitted))
names(singlet_scores) <- NULL

# Fill the object as run.RCTD(doublet_mode = "doublet") would (no score_mat from rctd-py)
RCTD@internal_vars$sigma <- run_info$sigma / 100
RCTD@config$RCTDmode     <- "doublet"
RCTD@results <- list(
  results_df      = results_df,
  weights         = weights,
  weights_doublet = weights_doublet,
  score_mat       = NULL,
  singlet_scores  = singlet_scores
)

# From here on, the outputs of rctd.R
weights_all <- as.data.frame(as.matrix(RCTD@results$weights))[xe_cells, ]
weights_all$cell_id <- xe_cells
weights_all <- weights_all %>% select(cell_id, everything())

annotations.df <- RCTD@results$results_df[xe_cells, ]
annotations.df <- annotations.df %>% mutate(second_type = scond_type)
annotations.df$scond_type <- NULL
labels         <- annotations.df %>% select(first_type, second_type)
labels$cell_id <- xe_cells
labels         <- labels %>% select(cell_id, everything())

# Second type of highly confident cells, computed as in rctd.R (weights data.frame incl. cell_id)
weights_bool <- weights_all > 0.01
n_candidates <- rowSums(weights_bool)
high_confidence_cells <- n_candidates < 2
labels <- labels %>%
  mutate(second_type_updated = if_else(high_confidence_cells, NA_character_, second_type))

saveRDS(RCTD, snakemake@output[["rds_output"]])
write_parquet(labels, snakemake@output[["labels"]])
write_parquet(weights_all, snakemake@output[["scores"]])

results_out <- RCTD@results$results_df
results_out$cell_id <- rownames(results_out)
results_out <- results_out %>% select(cell_id, everything())

write_parquet(results_out, snakemake@output[["out_res_df"]])
write_parquet(as.data.frame(as.matrix(RCTD@results$weights)), snakemake@output[["out_w"]])
write_parquet(as.data.frame(as.matrix(RCTD@results$weights_doublet)), snakemake@output[["out_wd"]])
# Singlet scores in the long format of convert_singlet_scores_to_long_df (cell_index = row of results_df)
write_parquet(
  data.frame(cell_index = match(py_sc$cell_id, fitted), cell_type = py_sc$cell_type, score = py_sc$score),
  snakemake@output[["out_sc"]]
)
# rctd-py returns no pair score matrices: an empty table with rctd.R's columns
write_parquet(data.frame(row = character(), column = character(), value = numeric(), cell_index = integer()), snakemake@output[["out_sm"]])
