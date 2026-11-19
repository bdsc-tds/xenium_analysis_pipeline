log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Identify marker genes per cell type using scran::scoreMarkers, and aggregate unnormalised counts per cell type (pseudobulk)

library(Seurat)
library(scran)
library(BiocParallel)
library(Matrix)
library(arrow)
library(dplyr)

###### snakemake params  ######

xe_assay           <- snakemake@params[["xe_assay"]] # assay with log-normalised expression used for marker detection
counts_assay       <- snakemake@params[["counts_assay"]] # assay with unnormalised counts used for the pseudobulk table (corrected counts if run after count correction)
control_gene_pat   <- snakemake@params[["control_gene_pat"]]
min_cells_per_type <- snakemake@params[["min_cells_per_type"]] %||% 25
rctd_singlets_only <- snakemake@params[["rctd_singlets_only"]] %||% TRUE
threads            <- snakemake@threads

rctd_results_path  <- snakemake@input[["rctd_results"]] # empty unless the annotation comes from RCTD

###### end of snakemake params  ######

# Columns of the marker table besides the effect sizes, in the same order as returned by scran::scoreMarkers
marker_stats_cols <- c(
  "self.average", "other.average", "self.detected", "other.detected",
  as.vector(t(outer(
    c("logFC.cohen", "AUC", "logFC.detected"),
    c("mean", "min", "median", "max", "rank"),
    function(effect, summary) paste(summary, effect, sep = ".")
  )))
)

###### Load data  ######

message("Loading data...")
xe <- readRDS(snakemake@input[["xe"]])
labels_df <- as.data.frame(read_parquet(snakemake@input[["labels"]]))
message("Data loaded.")

# `Assays` is masked by SummarizedExperiment, hence the explicit namespace
for (assay in unique(c(xe_assay, counts_assay))) {
  if (!assay %in% as.character(Seurat::Assays(xe))) stop("Error: Assay '", assay, "' is not found in the xenium object.")
}

# RCTD stores the label in `first_type`, the other annotation methods in `label`
label_col <- intersect(c("label", "first_type"), colnames(labels_df))[1]
if (is.na(label_col)) stop("Error: Neither 'label' nor 'first_type' is found in the labels table.")
if (!"cell_id" %in% colnames(labels_df)) stop("Error: 'cell_id' is not found in the labels table.")

labels_df <- data.frame(
  cell_id = as.character(labels_df$cell_id),
  cell_type = as.character(labels_df[[label_col]])
)
message(nrow(labels_df), " cells in the labels table, ", ncol(xe), " cells in the xenium object.")

###### Select cells  ######

labels_df <- labels_df %>%
  filter(!is.na(cell_type), cell_id %in% colnames(xe))
message(nrow(labels_df), " labelled cells found in the xenium object.")

if (length(rctd_results_path) > 0 && rctd_singlets_only) {
  rctd_df <- as.data.frame(read_parquet(rctd_results_path))
  if (!all(c("cell_id", "spot_class") %in% colnames(rctd_df))) {
    stop("Error: 'cell_id' and 'spot_class' are required in the RCTD results table to select singlets.")
  }

  singlets <- rctd_df$cell_id[!is.na(rctd_df$spot_class) & rctd_df$spot_class == "singlet"]
  labels_df <- labels_df %>% filter(cell_id %in% singlets)
  message(nrow(labels_df), " cells left after keeping RCTD singlets only.")
}

###### Extract expression  ######

expr   <- as(GetAssayData(xe, assay = xe_assay, layer = "data"), "CsparseMatrix")
counts <- as(GetAssayData(xe, assay = counts_assay, layer = "counts"), "CsparseMatrix")
rm(xe)
invisible(gc())

# Control probes are not expected in these assays, but make sure they are excluded
# Seurat replaces underscores in feature names with dashes, hence both forms are checked
if (!is.null(control_gene_pat) && nzchar(control_gene_pat)) {
  is_control <- function(genes) grepl(control_gene_pat, genes) | grepl(control_gene_pat, gsub("-", "_", genes))

  expr   <- expr[!is_control(rownames(expr)), , drop = FALSE]
  counts <- counts[!is_control(rownames(counts)), , drop = FALSE]
}

expr   <- expr[, labels_df$cell_id, drop = FALSE]
counts <- counts[, labels_df$cell_id, drop = FALSE]

type_sizes <- table(labels_df$cell_type)
message("Number of cells per cell type:")
print(type_sizes)

###### Pseudobulk  ######
# All the labelled cell types are kept regardless of their sizes, so that filtering can be done after gathering samples

if (nrow(labels_df) > 0) {
  message("Aggregating unnormalised counts per cell type...")
  type_indicator <- Matrix::t(Matrix::fac2sparse(factor(labels_df$cell_type))) # cells x cell types

  detected <- counts
  detected@x <- as.numeric(detected@x > 0)

  sums       <- as.matrix(counts %*% type_indicator) # genes x cell types
  n_detected <- as.matrix(detected %*% type_indicator)

  pseudobulk <- data.frame(
    cell_type = rep(colnames(sums), each = nrow(sums)),
    gene = rep(rownames(sums), times = ncol(sums)),
    counts = as.vector(sums),
    n_cells = rep(as.integer(type_sizes[colnames(sums)]), each = nrow(sums)),
    n_detected = as.integer(round(as.vector(n_detected)))
  )
  rm(detected, type_indicator)
} else {
  message("Warning: No labelled cells available. An empty pseudobulk table will be saved.")
  pseudobulk <- data.frame(
    cell_type = character(0),
    gene = character(0),
    counts = numeric(0),
    n_cells = integer(0),
    n_detected = integer(0)
  )
}
rm(counts)
invisible(gc())

###### Marker genes  ######

kept_types    <- names(type_sizes)[type_sizes >= min_cells_per_type]
dropped_types <- setdiff(names(type_sizes), kept_types)
if (length(dropped_types) > 0) {
  message(
    "Cell types with fewer than ", min_cells_per_type, " cells are excluded from marker detection: ",
    paste0(dropped_types, " (", type_sizes[dropped_types], ")", collapse = ", ")
  )
}

if (length(kept_types) >= 2) {
  is_kept <- labels_df$cell_type %in% kept_types

  message("Running scoreMarkers on ", sum(is_kept), " cells of ", length(kept_types), " cell types with ", threads, " threads...")
  scores <- scoreMarkers(
    expr[, is_kept, drop = FALSE],
    groups = labels_df$cell_type[is_kept],
    BPPARAM = if (threads > 1) MulticoreParam(workers = threads) else SerialParam()
  )
  message("Done running scoreMarkers.")

  markers <- lapply(names(scores), function(cell_type) {
    stats <- as.data.frame(scores[[cell_type]])
    data.frame(
      cell_type = cell_type,
      gene = rownames(stats),
      n_cells = as.integer(type_sizes[[cell_type]]),
      stats,
      check.names = FALSE,
      row.names = NULL
    )
  }) %>% bind_rows()
} else {
  message("Warning: Fewer than two cell types with at least ", min_cells_per_type, " cells. An empty marker table will be saved.")
  markers <- data.frame(
    cell_type = character(0),
    gene = character(0),
    n_cells = integer(0)
  )
  for (col in marker_stats_cols) markers[[col]] <- numeric(0)
}

###### Save results  ######

write_parquet(markers, snakemake@output[["markers"]])
write_parquet(pseudobulk, snakemake@output[["pseudobulk"]])
message("Results saved.")
