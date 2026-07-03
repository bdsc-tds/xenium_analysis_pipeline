# =============================================================================
# seurat_to_h5_parquet.R  (Snakemake `script:` entry point)
#
# Convert a Seurat scRNA-seq object to:
#   - <h5_path>       10x-compatible HDF5 count matrix
#   - <parquet_path>  cell metadata (+ optional embeddings)
#
# Expects to be invoked via a Snakemake `script:` directive from a rule
# shaped like `runConvertSeuratReferenceToH5` (see seurat2h5.smk), which
# injects the `snakemake` S4 object with:
#
#   input$rds               path to the Seurat .rds object
#   output$h5                target .h5 path
#   output$parquet           target .parquet path
#   params$assay             Seurat assay to export (default "RNA")
#   params$layer             Seurat layer/slot to export (default "counts")
#   params$genome            genome tag written into the 10x HDF5 (default "unknown")
#   log[[1]]                 log file path
# =============================================================================

## ---- log setup -------------------------------------------------------
log_file <- snakemake@log[[1]]
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(arrow)
})

# =============================================================================
# Main conversion function
# =============================================================================
seurat_to_h5_parquet <- function(seurat_obj,
                                 h5_path,
                                 parquet_path,
                                 assay            = "RNA",
                                 layer            = "counts",
                                 version          = "3",
                                 genome           = "unknown",
                                 extra_reductions = TRUE,
                                 overwrite        = TRUE) {
  
  if (!exists("write10xCounts", mode = "function")) {
    stop(
      "write10xCounts() not found.\n",
      "Source the bdsc-tds readwrite.R before calling this function."
    )
  }
  
  stopifnot(
    "seurat_obj must be a Seurat object" = inherits(seurat_obj, "Seurat"),
    "h5_path must be a character string" = is.character(h5_path),
    "parquet_path must be a character string" = is.character(parquet_path),
    "assay not found in Seurat object"   = assay %in% Assays(seurat_obj)
  )
  
  for (out_path in c(h5_path, parquet_path)) {
    out_dir <- dirname(out_path)
    if (!dir.exists(out_dir) && nzchar(out_dir) && out_dir != ".") {
      message(sprintf("Creating output directory: %s", out_dir))
      dir.create(out_dir, recursive = TRUE)
    }
  }
  
  DefaultAssay(seurat_obj) <- assay
  
  # ------------------------------------------------------------------
  # 1. Extract count matrix (Seurat v4 and v5 compatible)
  # ------------------------------------------------------------------
  seurat_v5 <- utils::packageVersion("Seurat") >= "5.0.0"
  
  if (seurat_v5) {
    available_layers <- Layers(seurat_obj[[assay]])
    if (!layer %in% available_layers) {
      stop(sprintf(
        "Layer '%s' not found in assay '%s'.\nAvailable: %s",
        layer, assay, paste(available_layers, collapse = ", ")
      ))
    }
    counts_mat <- LayerData(seurat_obj, assay = assay, layer = layer)
  } else {
    counts_mat <- GetAssayData(seurat_obj, assay = assay, slot = layer)
  }
  
  message(sprintf(
    "[seurat_to_h5_parquet] assay='%s'  layer='%s'  (%d genes x %d cells)",
    assay, layer, nrow(counts_mat), ncol(counts_mat)
  ))
  
  # ------------------------------------------------------------------
  # 2. Write counts -> HDF5 (delegated to readwrite.R)
  # ------------------------------------------------------------------
  message(sprintf("  Writing HDF5    -> %s", h5_path))
  
  if (overwrite && file.exists(h5_path)) file.remove(h5_path)
  
  write10xCounts(
    path      = h5_path,
    x         = counts_mat,
    type      = "HDF5",
    genome    = genome,
    version   = version,
    overwrite = FALSE    # already handled above
  )
  
  # ------------------------------------------------------------------
  # 3. Build metadata data.frame
  # ------------------------------------------------------------------
  meta <- seurat_obj@meta.data
  
  if (extra_reductions && length(seurat_obj@reductions) > 0) {
    for (red_name in names(seurat_obj@reductions)) {
      emb    <- Embeddings(seurat_obj, reduction = red_name)
      n_keep <- min(ncol(emb), 10L)
      emb    <- emb[, seq_len(n_keep), drop = FALSE]
      colnames(emb) <- paste0(red_name, "_", seq_len(n_keep))
      shared <- intersect(rownames(meta), rownames(emb))
      if (length(shared)) {
        meta[shared, colnames(emb)] <- emb[shared, , drop = FALSE]
        message(sprintf(
          "  Appended reduction '%s' (%d dims) to metadata.", red_name, n_keep
        ))
      }
    }
  }
  
  # ------------------------------------------------------------------
  # 4. Write metadata -> Parquet
  # ------------------------------------------------------------------
  message(sprintf("  Writing Parquet -> %s", parquet_path))
  
  df <- as.data.frame(meta)
  df[["cell_id"]] <- rownames(df)
  df <- df[, c("cell_id", setdiff(colnames(df), "cell_id")), drop = FALSE]
  rownames(df) <- NULL
  
  for (col in colnames(df)) {
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
  }
  
  arrow::write_parquet(df, parquet_path)
  
  # ------------------------------------------------------------------
  # 5. Summary
  # ------------------------------------------------------------------
  message("\n============== Conversion complete ==============")
  message(sprintf("  HDF5    : %s", h5_path))
  message(sprintf("  Parquet : %s", parquet_path))
  message(sprintf("  Cells   : %d", ncol(counts_mat)))
  message(sprintf("  Genes   : %d", nrow(counts_mat)))
  message(sprintf("  Metadata columns : %d", ncol(meta)))
  message("=================================================\n")
  
  invisible(list(h5 = h5_path, parquet = parquet_path))
}

# =============================================================================
# Snakemake entry point
# =============================================================================
input_path     <- snakemake@input[["rds"]]
h5_path        <- snakemake@output[["h5"]]
parquet_path   <- snakemake@output[["parquet"]]
assay          <- snakemake@params[["assay"]]
layer          <- snakemake@params[["layer"]]
genome         <- if (!is.null(snakemake@params[["genome"]])) snakemake@params[["genome"]] else "unknown"

snakemake@source("../../scripts/utils/readwrite.R")

message(sprintf("[snakemake] Loading Seurat object:    %s", input_path))
seurat_obj <- readRDS(input_path)

seurat_to_h5_parquet(
  seurat_obj   = seurat_obj,
  h5_path      = h5_path,
  parquet_path = parquet_path,
  assay        = assay,
  layer        = layer,
  genome       = genome
)

sink(type = "message")
sink(type = "output")

