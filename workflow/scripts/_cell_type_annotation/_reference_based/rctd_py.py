"""
RCTD in doublet mode with rctd-py, reproducing the reference filtering and settings of rctd.R
(generate_reference_obj, generate_class_df and create.RCTD in annotation_utils.R / rctd.R).
The R rule assembleReferenceBasedRCTDPy turns the outputs into the files rctd.R writes.
"""

import argparse
import importlib.metadata
import json

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from rctd import RCTDConfig, Reference
from rctd._doublet import run_doublet_mode
from rctd._irwls import solve_irwls_batch
from rctd._likelihood import calc_log_likelihood_batch
from rctd._rctd import RCTD
from rctd._types import SPOT_CLASS_NAMES


def read_10x_h5(path: str) -> anndata.AnnData:
    """Read a 10x v3 HDF5 count matrix (genes x cells, CSC) into a cells x genes AnnData."""
    with h5py.File(path, "r") as f:
        # Build the CSC matrix stored by write10xCounts
        g = f["matrix"]
        n_genes, n_cells = g["shape"][()]
        mat = sp.csc_matrix(
            (g["data"][()], g["indices"][()], g["indptr"][()]),
            shape=(n_genes, n_cells),
        )
        barcodes = [b.decode() for b in g["barcodes"][()]]
        genes = [x.decode() for x in g["features"]["name"][()]]

    # Cells x genes, as AnnData expects
    adata = anndata.AnnData(X=mat.T.tocsr())
    adata.obs_names = barcodes
    adata.var_names = genes
    return adata


def filter_reference(
    ref: anndata.AnnData,
    meta: pd.DataFrame,
    annotation_level: str,
    query_genes: list[str],
    donor: str,
    min_umi: int,
    max_umi: int,
    cell_min_instance: int,
) -> anndata.AnnData:
    """generate_reference_obj (annotation_utils.R): donor-specific types, UMI range, rare types."""
    # Align metadata to the reference cells
    meta = meta.set_index("cell_id").loc[ref.obs_names]
    if annotation_level not in meta.columns:
        raise ValueError(f"Annotation level {annotation_level} is missing in the reference metadata.")
    ref.obs = meta.copy()
    ref = ref[~ref.obs[annotation_level].isna()].copy()

    # Donor-specific cell types are kept only for the query's donor (labels matched as a regex, as in R)
    if "donor" in ref.obs.columns:
        ref_donors = pd.unique(ref.obs["donor"].astype(str))
        if donor in ref_donors:
            pattern = "|".join(ref_donors)
            is_donor_specific = ref.obs[annotation_level].astype(str).str.contains(pattern, regex=True)
            keep = (~is_donor_specific) | (ref.obs["donor"].astype(str) == donor)
            print(f"Donor {donor} is in the reference: removed {int((~keep).sum())} cells of other donors' cell types.")
            ref = ref[keep.values].copy()

    # UMI counts on the genes shared with the query, kept strictly between min_umi and max_umi
    features = [g for g in ref.var_names if g in set(query_genes)]
    n_count = np.asarray(ref[:, features].X.sum(axis=1)).ravel()
    ref = ref[(n_count > min_umi) & (n_count < max_umi)].copy()

    # Cell types with fewer than cell_min_instance cells are removed
    type_counts = ref.obs[annotation_level].astype(str).value_counts()
    keep_types = type_counts.index[type_counts >= cell_min_instance]
    dropped = sorted(set(type_counts.index) - set(keep_types))
    if dropped:
        print(f"Cell types {','.join(dropped)} were filtered out from the reference due to low abundance.")
    ref = ref[ref.obs[annotation_level].astype(str).isin(keep_types).values].copy()
    return ref


def make_class_df(obs: pd.DataFrame, annotation_level: str, class_level: str | None) -> dict[str, str] | None:
    """generate_class_df (annotation_utils.R): most frequent class of each cell type."""
    if not class_level or class_level == annotation_level:
        return None
    if class_level not in obs.columns:
        print(f"class_level '{class_level}' is not in the reference metadata; running with class_df = None.")
        return None

    # Count (type, class) pairs; the most frequent class wins, ties broken by class order as dplyr does
    counts = (
        obs[[annotation_level, class_level]]
        .astype(str)
        .value_counts()
        .rename("n")
        .reset_index()
        .sort_values([annotation_level, "n", class_level], ascending=[True, False, True])
    )
    best = counts.drop_duplicates(annotation_level, keep="first")
    return dict(zip(best[annotation_level], best[class_level]))


def singlet_scores(rctd: RCTD, full_weights: np.ndarray, batch_size: int, device: str) -> pd.DataFrame:
    """
    Singlet score of every candidate type of every cell (spacexr's results$singlet_scores), as in
    step 4 of rctd-py's run_doublet_mode, which does not return them. SPLIT's post-processing needs them.
    Returns the long table of rctd.R (cell_index 1-based, cell_type, score).
    """
    n_cells, n_types = full_weights.shape

    # Candidate types: full-mode weight above 0.01, with spacexr's fallbacks for 0 or 1 candidate
    singles = []
    for n in range(n_cells):
        cands = np.where(full_weights[n] > 0.01)[0].tolist()
        if len(cands) == 0:
            cands = list(range(min(3, n_types)))
        elif len(cands) == 1:
            cands.append(1 if cands[0] == 0 else 0)
        singles.extend((n, t) for t in cands)
    singles = np.asarray(singles, dtype=np.int64)

    # Data and model tensors on the device
    P = torch.tensor(rctd.norm_profiles, device=device)
    Q = torch.tensor(rctd.q_mat, device=device)
    SQ = torch.tensor(rctd.sq_mat, device=device)
    X = torch.tensor(rctd.x_vals, device=device)
    Y = torch.tensor(rctd.counts, device=device)
    nUMI = torch.tensor(rctd.nUMI, device=device)

    # Fit each (cell, type) as a singlet and score it, in batches
    scores = []
    for start in range(0, len(singles), batch_size):
        sg = singles[start : start + batch_size]
        pix = torch.as_tensor(sg[:, 0], device=device)
        typ = torch.as_tensor(sg[:, 1], device=device)
        S = nUMI[pix][:, None, None] * P[:, typ].T.unsqueeze(-1)
        w, _ = solve_irwls_batch(
            S_batch=S,
            Y_batch=Y[pix],
            nUMI_batch=nUMI[pix],
            Q_mat=Q,
            SQ_mat=SQ,
            x_vals=X,
            max_iter=25,
            min_change=0.001,
            constrain=False,
            bulk_mode=False,
        )
        expected = torch.clamp(torch.sum(S * w[:, None, :], dim=-1), min=1e-4)
        scores.append(calc_log_likelihood_batch(Y[pix], expected, Q, SQ, X))

    # Long table in cell order, candidates in type order
    out = pd.DataFrame(
        {
            "cell_index": singles[:, 0] + 1,
            "cell_type": np.asarray(rctd.reference.cell_type_names)[singles[:, 1]],
            "score": torch.cat(scores).cpu().numpy().astype(float),
        }
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref_counts", required=True)
    parser.add_argument("--ref_metadata", required=True)
    parser.add_argument("--query_counts", required=True)
    parser.add_argument("--query_metadata", required=True)
    parser.add_argument("--annotation_level", required=True)
    parser.add_argument("--class_level", default="")
    parser.add_argument("--donor", required=True)
    parser.add_argument("--ref_min_umi", type=int, required=True)
    parser.add_argument("--ref_max_umi", type=int, required=True)
    parser.add_argument("--xe_min_umi", type=int, required=True)
    parser.add_argument("--xe_min_counts", type=int, required=True)
    parser.add_argument("--cell_min_instance", type=int, required=True)
    parser.add_argument("--umi_min_sigma", type=int, required=True)
    parser.add_argument("--batch_size", type=int, default=10000)
    parser.add_argument("--dtype", default="float64", choices=["float64", "float32"])
    parser.add_argument("--is_int_counts", action="store_true")
    parser.add_argument("--out_results_df", required=True)
    parser.add_argument("--out_weights", required=True)
    parser.add_argument("--out_weights_doublet", required=True)
    parser.add_argument("--out_singlet_scores", required=True)
    parser.add_argument("--out_reference_cells", required=True)
    parser.add_argument("--out_class_df", required=True)
    parser.add_argument("--out_run_info", required=True)
    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"rctd-py {importlib.metadata.version('rctd-py')}, torch {torch.__version__}, device {device}")

    # Load the query (all genes, as SpatialRNA in rctd.R)
    query = read_10x_h5(args.query_counts)
    query_meta = pd.read_parquet(args.query_metadata).set_index("cell_id")
    if not query_meta.index.equals(query.obs_names):
        raise ValueError("Cell ids of the query counts and metadata differ.")
    if not args.is_int_counts:
        print("Query counts are not integers (proseg expected counts); rctd-py uses them as given.")
    print(f"Query: {query.n_obs} cells x {query.n_vars} genes")

    # Load and filter the reference as rctd.R does
    ref = read_10x_h5(args.ref_counts)
    ref_meta = pd.read_parquet(args.ref_metadata)
    ref = filter_reference(
        ref,
        ref_meta,
        args.annotation_level,
        list(query.var_names),
        args.donor,
        args.ref_min_umi,
        args.ref_max_umi,
        args.cell_min_instance,
    )
    class_df = make_class_df(ref.obs, args.annotation_level, args.class_level)

    # Reference profiles on the genes shared with the query (Reference() in rctd.R)
    common_genes = [g for g in query.var_names if g in set(ref.var_names)]
    ref = ref[:, common_genes].copy()
    ref.obs["cell_type"] = ref.obs[args.annotation_level].astype(str).values
    reference = Reference(
        ref,
        cell_type_col="cell_type",
        cell_min=args.cell_min_instance,
        min_UMI=args.ref_min_umi,
    )
    print(f"Reference: {ref.n_obs} cells, {reference.n_types} cell types, {len(common_genes)} shared genes")

    # create.RCTD settings of rctd.R; no torch.compile (slow first compile with this torch version)
    config = RCTDConfig(
        UMI_min=args.xe_min_umi,
        counts_MIN=args.xe_min_counts,
        UMI_min_sigma=args.umi_min_sigma,
        class_df=class_df,
        dtype=args.dtype,
        device=device,
        compile=False,
    )

    # Fit the platform effect and sigma, then doublet mode (run.RCTD(doublet_mode = "doublet"))
    rctd = RCTD(query, reference, config)
    rctd.fit_platform_effects()
    res = run_doublet_mode(
        spatial_counts=rctd.counts,
        spatial_numi=rctd.nUMI,
        norm_profiles=rctd.norm_profiles,
        cell_type_names=rctd.reference.cell_type_names,
        q_mat=rctd.q_mat,
        sq_mat=rctd.sq_mat,
        x_vals=rctd.x_vals,
        config=rctd.config,
        batch_size=args.batch_size,
        device=device,
        pixel_mask=rctd._pixel_mask,
    )

    # Cells kept by the UMI and counts_MIN filters, row order of every result array
    cell_ids = np.asarray(query.obs_names)[rctd._pixel_mask]
    types = list(res.cell_type_names)
    print(f"Fitted {len(cell_ids)} of {query.n_obs} cells, sigma {rctd.sigma}")

    # results_df as spacexr names it (scond_type)
    results_df = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "spot_class": [SPOT_CLASS_NAMES[i] for i in res.spot_class],
            "first_type": [types[i] for i in res.first_type],
            "scond_type": [types[i] for i in res.second_type],
            "first_class": np.asarray(res.first_class, dtype=bool),
            "second_class": np.asarray(res.second_class, dtype=bool),
            "min_score": np.asarray(res.min_score, dtype=float),
            "singlet_score": np.asarray(res.singlet_score, dtype=float),
        }
    )
    results_df.to_parquet(args.out_results_df, index=False)

    # Full-mode weights and doublet weights, one row per fitted cell
    weights = pd.DataFrame(np.asarray(res.weights, dtype=float), columns=types)
    weights.insert(0, "cell_id", cell_ids)
    weights.to_parquet(args.out_weights, index=False)
    weights_doublet = pd.DataFrame(np.asarray(res.weights_doublet, dtype=float), columns=["first_type", "second_type"])
    weights_doublet.insert(0, "cell_id", cell_ids)
    weights_doublet.to_parquet(args.out_weights_doublet, index=False)

    # Per-type singlet scores of the candidate types (results$singlet_scores in spacexr)
    singlet_scores(rctd, np.asarray(res.weights), args.batch_size, device).to_parquet(
        args.out_singlet_scores, index=False
    )

    # Reference cells and labels after filtering, so R builds the same reference object
    pd.DataFrame({"cell_id": ref.obs_names, "cell_type": ref.obs["cell_type"].values}).to_parquet(
        args.out_reference_cells, index=False
    )
    pd.DataFrame(
        {"cell_type": list(class_df.keys()) if class_df else [], "class": list(class_df.values()) if class_df else []}
    ).to_parquet(args.out_class_df, index=False)

    # Run details for the R side and for the record
    run_info = {
        "rctd_py_version": importlib.metadata.version("rctd-py"),
        "torch_version": torch.__version__,
        "device": device,
        "sigma": int(rctd.sigma),
        "cell_type_names": types,
        "n_query_cells": int(query.n_obs),
        "n_fitted_cells": int(len(cell_ids)),
        "n_reference_cells": int(ref.n_obs),
        "batch_size": args.batch_size,
        "dtype": args.dtype,
        "class_df": class_df is not None,
    }
    with open(args.out_run_info, "w") as f:
        json.dump(run_info, f, indent=2)


if __name__ == "__main__":
    main()
