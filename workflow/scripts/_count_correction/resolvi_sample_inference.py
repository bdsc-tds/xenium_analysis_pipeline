import dask

dask.config.set({"dataframe.query-planning": False})

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import argparse
import os
import sys

from .._joint_scanpy_analysis import preprocessing
from ..utils import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument("--path", type=str, 
        required=True, help="Path to the xenium sample file.",)
    parser.add_argument("--dir_resolvi_model", type=str, 
        required=True, help="directory with saved RESOLVI model weights",)
    parser.add_argument(
        "--out_file_resolvi_corrected_counts",
        type=str,
        help="Path to resolvi corrected counts parquet file.",
    )
    parser.add_argument(
        "--out_file_resolvi_proportions",
        type=str,
        help="Path to resolvi proportions parquet file.",
    )
    parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
    parser.add_argument(
        "--num_samples",
        type=int,
        help="Number of samples for RESOLVI generative model.",
    )
    parser.add_argument("--batch_size", type=int, default=1000, help="batch size parameter")
    parser.add_argument("--cell_type_labels", type=str, help="optional cell_type_labels for semi-supervised mode")

    ret = parser.parse_args()
    if not os.path.isdir(ret.path):
        raise RuntimeError(f"Error! Input directory to data does not exist: {ret.path}")
    if not os.path.isdir(ret.dir_resolvi_model):
        raise RuntimeError(f"Error! Input directory to trained model does not exist: {ret.dir_resolvi_model}")

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # read counts
    adata = readwrite.read_xenium_sample(
        args.path,
        cells_as_circles=False,
        cells_boundaries=False,
        cells_boundaries_layers=False,
        nucleus_boundaries=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
        aligned_images=False,
        anndata=True,
    )

    # need to round proseg expected counts for resolVI to run
    # no need for if statement, doesn't change anything to other segmentation methods
    adata.X.data = adata.X.data.astype(np.float32).round()
    adata.obs_names = adata.obs_names.astype(str)

    if args.cell_type_labels is not None:
        labels_key = "labels_key"
        semisupervised = True
        adata.obs[labels_key] = pd.read_parquet(args.cell_type_labels).iloc[:, 0]
    else:
        labels_key = None
        semisupervised = False

    # preprocess (QC filters only)
    # resolvi requires at least 5 counts in each cell
    preprocessing.preprocess(
        adata,
        normalize=False,
        log1p=False,
        scale="none",
        pca=False,
        umap=False,
        save_raw=False,
        min_counts=args.min_counts,
        min_genes=args.min_features,
        max_counts=args.max_counts,
        max_genes=args.max_features,
        min_cells=args.min_cells,
        backend="cpu",
    )

    resolvi = scvi.external.RESOLVI.load(args.dir_resolvi_model, adata=adata)

    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={"post_sample_q50": np.median},
        num_samples=args.num_samples,
        batch_size=args.batch_size,
        summary_frequency=100,
    )
    samples_corr = pd.DataFrame(samples_corr).T

    samples = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=["mixture_proportions"],
        summary_fun={"post_sample_means": np.mean},
        num_samples=args.num_samples,
        batch_size=args.batch_size,
        summary_frequency=100,
    )
    samples_proportions = pd.DataFrame(samples).T

    ### save
    samples_corr = pd.DataFrame(
        samples_corr.loc["post_sample_q50", "px_rate"],
        index=adata.obs_names,
        columns=adata.var_names,
    )
    samples_proportions = pd.DataFrame(
        samples_proportions.loc["post_sample_means", "mixture_proportions"],
        index=adata.obs_names,
        columns=["true_proportion", "diffusion_proportion", "background_proportion"],
    )

    adata_out = ad.AnnData(samples_corr)
    readwrite.write_10X_h5(adata_out, args.out_file_resolvi_corrected_counts)
    samples_proportions.to_parquet(args.out_file_resolvi_proportions)
    resolvi.save(args.out_dir_resolvi_model)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
