import dask

dask.config.set({"dataframe.query-planning": False})

import spatialdata_io
import numpy as np
import pandas as pd
import scvi
import sys
import argparse
import anndata as ad
from pathlib import Path

sys.path.append("workflow/scripts/")
import preprocessing
import readwrite

# params
parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")
parser.add_argument("--path", type=Path, help="Path to the xenium donor file.")
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
parser.add_argument(
    "--min_features", type=int, help="QC parameter from pipeline config"
)
parser.add_argument(
    "--max_counts", type=float, help="QC parameter from pipeline config"
)
parser.add_argument(
    "--max_features", type=float, help="QC parameter from pipeline config"
)
parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
parser.add_argument(
    "--num_samples", type=int, help="Number of samples for RESOLVI generative model."
)
args = parser.parse_args()

# Access the arguments
path = args.path
out_file_resolvi_corrected_counts = args.out_file_resolvi_corrected_counts
out_file_resolvi_proportions = args.out_file_resolvi_proportions
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells
num_samples = args.num_samples

# read counts
adata = spatialdata_io.xenium(
    path,
    cells_as_circles=False,
    cells_boundaries=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
)["table"]


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
    min_counts=min_counts,
    min_genes=min_features,
    max_counts=max_counts,
    max_genes=max_features,
    min_cells=min_cells,
    backend="cpu",
)


scvi.external.RESOLVI.setup_anndata(
    adata, labels_key=None, layer=None, prepare_data_kwargs={"spatial_rep": "spatial"}
)
resolvi = scvi.external.RESOLVI(adata, semisupervised=False)
resolvi.train(max_epochs=50)

samples_corr = resolvi.sample_posterior(
    model=resolvi.module.model_corrected,
    return_sites=["px_rate"],
    summary_fun={"post_donor_q50": np.median},
    num_samples=num_samples,
    summary_frequency=100,
)
samples_corr = pd.DataFrame(samples_corr).T

samples = resolvi.sample_posterior(
    model=resolvi.module.model_residuals,
    return_sites=["mixture_proportions"],
    summary_fun={"post_donor_means": np.mean},
    num_samples=num_samples,
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
readwrite.write_10X_h5(adata_out, out_file_resolvi_corrected_counts)
# samples_corr.to_parquet(out_file_resolvi_corrected)
samples_proportions.to_parquet(out_file_resolvi_proportions)
