import sys
import os
import argparse
import json

import pandas as pd
import scanpy as sc
import anndata as ad


import preprocessing


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")
    parser.add_argument(
        "-i",
        type=str,
        required=True,
        help="Path to the panel file.",
    )
    parser.add_argument(
        "-o",
        type=str,
        required=True,
        help="Path to the output file.",
    )
    parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument(
        "--n_comps",
        type=int,
        default=50,
        help="Number of components.",
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        default=50,
        help="Number of neighbors.",
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="cosine",
        help="Distance metric to use.",
    )
    parser.add_argument(
        "--min_dist",
        type=float,
        default=0.3,
        help="Minimum distance parameter.",
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        default=10,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--min_features",
        type=int,
        default=5,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--max_counts",
        type=float,
        default=float("inf"),
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--max_features",
        type=float,
        default=float("inf"),
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=5,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--use_gpu",
        action="store_true",
        help="whether to use gpu",
    )

    ret = parser.parse_args()
    if not os.path.isfile(ret.i):
        raise RuntimeError(f"Error! Input file does not exist: {ret.i}")

    os.makedirs(os.path.dirname(ret.o), exist_ok=True)

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    with open(
        args.i,
        mode="r",
        encoding="utf-8",
    ) as fh:
        files: dict[str, str] = json.load(fh)

    merged_ad: ad.AnnData = sc.concat(
        {k: sc.read_h5ad(v) for k, v in files.items()},
    )

    # preprocess
    preprocessing.preprocess(
        merged_ad,
        normalize=True,
        log1p=True,
        scale="none",
        n_comps=args.n_comps,
        metric=args.metric,
        min_dist=args.min_dist,
        n_neighbors=args.n_neighbors,
        backend="gpu" if args.use_gpu else "cpu",
        pca=True,
        umap=True,
        save_raw=False,
        min_counts=args.min_counts,
        min_genes=args.min_features,
        max_counts=args.max_counts,
        max_genes=args.max_features,
        min_cells=args.min_cells,
    )

    # save
    df_umap = pd.DataFrame(
        {
            "UMAP1": merged_ad.obsm["X_umap"][:, 0],
            "UMAP2": merged_ad.obsm["X_umap"][:, 1],
            "cell_id": merged_ad.obs["cell_id"],
            "sample_id": merged_ad.obs["sample_id"],
            "condition": merged_ad.obs["condition"],
            "gene_panel": merged_ad.obs["gene_panel"],
            "donor": merged_ad.obs["donor"],
            "sample": merged_ad.obs["sample"],
        },
        index=merged_ad.obs_names,
    )
    # df_umap[xenium_levels] = merged_ad.obs[xenium_levels]

    df_umap.to_parquet(args.o)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
