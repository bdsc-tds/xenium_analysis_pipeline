import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
import pandas as pd
import scvi
import argparse
import os
import sys

sys.path.append(
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
        ),
    ),
)
from _preprocessing import preprocessing
from utils import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(
        description="Run RESOLVI on a Xenium sample.",
    )
    parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument(
        "--path",
        type=str,
        required=True,
        help="Path to the xenium sample file.",
    )
    parser.add_argument(
        "--out_dir_resolvi_model",
        type=str,
        required=True,
        help="output directory with RESOLVI model weights",
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--min_features",
        type=int,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--max_counts",
        type=float,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--max_features",
        type=float,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        help="QC parameter from pipeline config",
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=50,
        help="Maximum number of epochs to train the model.",
    )
    parser.add_argument(
        "--cell_type_labels",
        type=str,
        help="optional cell_type_labels for semi-supervised mode",
    )
    parser.add_argument(
        "--mixture_k",
        type=int,
        default=50,
        help="mixture_k parameter for unsupervised RESOLVI",
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

    if args.cell_type_labels is not None:
        labels_key = "labels_key"
        semisupervised = True
        adata.obs[labels_key] = (
            pd.read_parquet(args.cell_type_labels)
            .set_index("cell_id")
            .iloc[:, 0]
            .astype("category")
        )
        adata = adata[adata.obs[labels_key].notna()]
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

    scvi.external.RESOLVI.setup_anndata(
        adata,
        labels_key=labels_key,
        layer=None,
        prepare_data_kwargs={"spatial_rep": "spatial"},
    )
    resolvi = scvi.external.RESOLVI(
        adata, mixture_k=args.mixture_k, semisupervised=semisupervised
    )
    resolvi.train(max_epochs=args.max_epochs)
    resolvi.save(args.out_dir_resolvi_model)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
