"""cell2location deconvolution of Xenium spatial transcriptomics data.

Reads Seurat RDS files via rpy2 for both query (Xenium) and reference (scRNA-seq),
trains a cell2location model, and writes labels/scores/abundance as parquet plus
the full AnnData as h5ad.

CLI mirrors the argparse pattern used by other Python scripts in this pipeline
(e.g. resolvi_sample_training.py). Run with:
    mamba run -n cell2location python3 cell2location.py --help
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run cell2location on a Xenium sample.")

    p.add_argument("--query", required=True,
                   help="Path to preprocessed_seurat.rds (Xenium Seurat object).")
    p.add_argument("--reference", required=True,
                   help="Path to reference scRNA-seq Seurat RDS.")
    p.add_argument("--annotation_id", required=True,
                   help="Full annotation_id wildcard (approach/ref_type/method/level/mode).")
    p.add_argument("--xe_assay", default="Xenium",
                   help="Assay name to pull counts from in the Xenium Seurat object.")
    p.add_argument("--ref_assay", default="RNA",
                   help="Assay name to pull counts from in the reference Seurat object.")

    p.add_argument("--out_h5ad", required=True, help="Output AnnData (h5ad).")
    p.add_argument("--out_labels", required=True, help="Output dominant cell type per cell (parquet).")
    p.add_argument("--out_scores", required=True, help="Output dominant cell type proportion per cell (parquet).")
    p.add_argument("--out_cell_abundance", required=True, help="Output full cell abundance matrix (parquet).")

    p.add_argument("--max_epochs_train", type=int, default=30000,
                   help="Training epochs for the NB regression reference model.")
    p.add_argument("--max_epochs_predict", type=int, default=1000,
                   help="Training epochs for the cell2location spatial model.")
    p.add_argument("--N_cells_per_location", type=float, default=1.0,
                   help="Expected number of cells per Xenium cell (typically 1 for single-cell resolution data).")
    p.add_argument("--detection_alpha", type=float, default=200.0,
                   help="Regularisation of the NegBinomial observation model.")

    p.add_argument("-l", dest="log", default=None, help="Path to log file.")

    return p.parse_args()


def _annotation_level(annotation_id: str) -> str:
    """Extract layer 3 (annotation_level) from the annotation_id path string."""
    parts = annotation_id.split("/")
    assert len(parts) == 5, f"Unexpected annotation_id format: {annotation_id}"
    return parts[3]


# ---------------------------------------------------------------------------
# Seurat → AnnData via rpy2
# ---------------------------------------------------------------------------

def _load_seurat_as_anndata(
    path: str,
    assay: str,
    label_col: str | None = None,
) -> ad.AnnData:
    """Read a Seurat RDS and return an AnnData with raw counts.

    The count matrix is extracted from the requested assay's 'counts' layer.
    If label_col is given, that metadata column is added to adata.obs['cell_type'].
    """
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    log = logging.getLogger(__name__)
    log.info("Loading %s (assay=%s) ...", path, assay)

    base = importr("base")
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(Matrix))")

    seurat_obj = base.readRDS(path)
    ro.globalenv["xe_obj"] = seurat_obj
    ro.globalenv["xe_assay"] = assay

    # Extract raw counts as a COO sparse matrix (genes × cells in R)
    ro.r("xe_mat <- GetAssayData(xe_obj, assay=xe_assay, layer='counts')")
    ro.r("xe_coo <- as(xe_mat, 'TsparseMatrix')")

    i_idx = np.asarray(ro.r("xe_coo@i"), dtype=np.int32)   # 0-based gene indices
    j_idx = np.asarray(ro.r("xe_coo@j"), dtype=np.int32)   # 0-based cell indices
    # Round to nearest integer: proseg and some other segmenters emit float counts
    vals  = np.round(np.asarray(ro.r("xe_coo@x"), dtype=np.float64)).astype(np.float32)
    nrow  = int(ro.r("nrow(xe_coo)")[0])
    ncol  = int(ro.r("ncol(xe_coo)")[0])

    # Build cells × genes CSR matrix
    counts = sp.coo_matrix((vals, (j_idx, i_idx)), shape=(ncol, nrow)).tocsr()

    gene_names = list(ro.r("rownames(xe_mat)"))
    cell_ids   = list(ro.r("colnames(xe_mat)"))

    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(index=cell_ids),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obs_names.name = "cell_id"
    adata.var_names.name = "gene"

    if label_col is not None:
        ro.globalenv["xe_label_col"] = label_col
        labels = list(ro.r("as.character(xe_obj@meta.data[[xe_label_col]])"))
        label_cell_ids = list(ro.r("rownames(xe_obj@meta.data)"))
        label_series = pd.Series(labels, index=label_cell_ids, name="cell_type")
        adata.obs["cell_type"] = label_series.reindex(adata.obs_names).values

    log.info("Loaded %d cells × %d genes.", adata.n_obs, adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    import cell2location
    import torch

    log = logging.getLogger(__name__)
    log.info("GPU available: %s", torch.cuda.is_available())

    annotation_level = _annotation_level(args.annotation_id)
    log.info("annotation_level = %s", annotation_level)

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    adata_xe = _load_seurat_as_anndata(args.query, assay=args.xe_assay)
    adata_ref = _load_seurat_as_anndata(
        args.reference, assay=args.ref_assay, label_col=annotation_level
    )

    # Drop reference cells with missing labels
    adata_ref = adata_ref[adata_ref.obs["cell_type"].notna()].copy()

    # ------------------------------------------------------------------
    # 2. Restrict to common genes
    # ------------------------------------------------------------------
    common_genes = adata_xe.var_names.intersection(adata_ref.var_names)
    log.info("Common genes: %d", len(common_genes))
    adata_xe  = adata_xe[:, common_genes].copy()
    adata_ref = adata_ref[:, common_genes].copy()

    # ------------------------------------------------------------------
    # 3. Train NB regression model on reference to learn gene signatures
    # ------------------------------------------------------------------
    log.info("Training RegressionModel (%d epochs) ...", args.max_epochs_train)
    cell2location.models.RegressionModel.setup_anndata(
        adata_ref, labels_key="cell_type"
    )
    reg_model = cell2location.models.RegressionModel(adata_ref)
    reg_model.train(
        max_epochs=args.max_epochs_train,
        accelerator="gpu",
    )

    adata_ref = reg_model.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": "gpu"},
    )
    # inf_aver: genes × cell_types (varm shape); Cell2location expects this orientation
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"].copy()
    else:
        inf_aver = adata_ref.varm["q05_cell_abundance_w_sf"].copy()
    log.info("Learned signatures for %d cell types.", inf_aver.shape[1])

    # ------------------------------------------------------------------
    # 4. Run cell2location spatial model on Xenium counts
    # ------------------------------------------------------------------
    log.info("Training Cell2location model (%d epochs) ...", args.max_epochs_predict)
    cell2location.models.Cell2location.setup_anndata(adata_xe)
    c2l_model = cell2location.models.Cell2location(
        adata_xe,
        cell_state_df=inf_aver,
        N_cells_per_location=args.N_cells_per_location,
        detection_alpha=args.detection_alpha,
    )
    c2l_model.train(
        max_epochs=args.max_epochs_predict,
        accelerator="gpu",
    )

    adata_xe = c2l_model.export_posterior(
        adata_xe,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": "gpu"},
    )

    # ------------------------------------------------------------------
    # 5. Extract results
    # ------------------------------------------------------------------
    abundance = pd.DataFrame(
        adata_xe.obsm["q05_cell_abundance_w_sf"],
        index=adata_xe.obs_names,
        columns=inf_aver.columns,
    )

    labels_df = (
        abundance.idxmax(axis=1, skipna=True)
        .rename("label")
        .reset_index()
    )
    labels_df.columns = ["cell_id", "label"]

    scores_df = (
        abundance.max(axis=1)
        .rename("score")
        .reset_index()
    )
    scores_df.columns = ["cell_id", "score"]

    abundance_df = abundance.reset_index().rename(columns={"cell_id": "cell_id"})

    # ------------------------------------------------------------------
    # 6. Save outputs
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(args.out_h5ad), exist_ok=True)

    log.info("Writing h5ad ...")
    adata_xe.write_h5ad(args.out_h5ad)

    log.info("Writing parquet outputs ...")
    labels_df.to_parquet(args.out_labels, index=False)
    scores_df.to_parquet(args.out_scores, index=False)
    abundance_df.to_parquet(args.out_cell_abundance, index=False)

    log.info("Done.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    _log_stream = None
    if args.log is not None:
        os.makedirs(os.path.dirname(args.log), exist_ok=True)
        _log_stream = open(args.log, "w", encoding="utf-8")
        sys.stdout = _log_stream
        sys.stderr = _log_stream

    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    try:
        main(args)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        if _log_stream is not None:
            _log_stream.close()
