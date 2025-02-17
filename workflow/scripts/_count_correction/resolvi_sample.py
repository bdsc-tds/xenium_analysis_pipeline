import os
import sys
import argparse

import numpy as np
import pandas as pd

import scvi
import anndata as ad

import dask

from .._joint_scanpy_analysis import preprocessing
from ..utils import readwrite


dask.config.set({"dataframe.query-planning": False})


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")
    parser.add_argument("--path", type=str, help="Path to the xenium donor file.")
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
    parser.add_argument(
        "--min_counts", type=int, help="QC parameter from pipeline config"
    )
    parser.add_argument(
        "--min_features", type=int, help="QC parameter from pipeline config"
    )
    parser.add_argument(
        "--max_counts", type=float, help="QC parameter from pipeline config"
    )
    parser.add_argument(
        "--max_features", type=float, help="QC parameter from pipeline config"
    )
    parser.add_argument(
        "--min_cells", type=int, help="QC parameter from pipeline config"
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=50,
        help="Maximum number of epochs to train the model.",
    )
    parser.add_argument(
        "--num_samples",
        type=int,
        help="Number of samples for RESOLVI generative model.",
    )

    ret = parser.parse_args()
    if not os.path.isdir(ret.path):
        raise RuntimeError(f"Error! Input directory does not exist: {ret.path}")
    if ret.max_epochs <= 0:
        ret.max_epochs = 50

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

    scvi.external.RESOLVI.setup_anndata(
        adata,
        labels_key=None,
        layer=None,
        prepare_data_kwargs={"spatial_rep": "spatial"},
    )
    resolvi = scvi.external.RESOLVI(adata, semisupervised=False)
    resolvi.train(max_epochs=args.max_epochs)

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

    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={"post_donor_q50": np.median},
        num_samples=args.num_samples,
        summary_frequency=100,
    )
    samples_corr = pd.DataFrame(samples_corr).T

    samples = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=["mixture_proportions"],
        summary_fun={"post_donor_means": np.mean},
        num_samples=args.num_samples,
        summary_frequency=100,
    )
    samples_proportions = pd.DataFrame(samples).T

    ### save
    samples_corr = pd.DataFrame(
        samples_corr.loc["post_donor_q50", "px_rate"],
        index=adata.obs_names,
        columns=adata.var_names,
    )
    samples_proportions = pd.DataFrame(
        samples_proportions.loc["post_donor_means", "mixture_proportions"],
        index=adata.obs_names,
        columns=["true_proportion", "diffusion_proportion", "background_proportion"],
    )

    adata_out = ad.AnnData(samples_corr)
    readwrite.write_10X_h5(adata_out, args.out_file_resolvi_corrected_counts)
    # samples_corr.to_parquet(out_file_resolvi_corrected)
    samples_proportions.to_parquet(args.out_file_resolvi_proportions)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
