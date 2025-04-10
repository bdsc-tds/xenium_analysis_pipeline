import dask

dask.config.set({"dataframe.query-planning": False})

from pathlib import Path
import argparse
import scipy
import pandas as pd
import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import preprocessing
import readwrite

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--sample", type=Path, help="Path to the sample file.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--xenium_processed_data_dir", type=Path, help="Path to the xenium processed data directories")
parser.add_argument("--normalisation_method", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--n_neighbors", type=int, help="Number of neighbors.")
parser.add_argument("--metric", type=str, help="Distance metric to use.")
parser.add_argument("--min_dist", type=float, help="Minimum distance parameter.")
parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
parser.add_argument("--genes", type=str, nargs="*", default=[], help="Restrict data to these genes for the UMAP.")

args = parser.parse_args()

# Access the arguments
xenium_processed_data_dir = args.xenium_processed_data_dir
sample = args.sample
out_file = args.out_file
normalisation_method = args.normalisation_method
layer = args.layer
n_comps = args.n_comps
n_neighbors = args.n_neighbors
metric = args.metric
min_dist = args.min_dist
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells
genes = args.genes


segmentation = sample.parents[3].stem
condition = sample.parents[2].stem
panel = sample.parents[1]
donor = sample.parents[0]

# read xenium sample
print("Reading sample")
print(donor.stem, sample.stem)

if segmentation == "proseg_expected":
    k = ("proseg", condition, panel.stem, donor.stem, sample.stem)
    name_sample = "/".join(k)
    sample_dir = xenium_processed_data_dir / f"{name_sample}/raw_results"
else:
    k = (segmentation.replace("proseg_mode", "proseg"), condition, panel.stem, donor.stem, sample.stem)
    name_sample = "/".join(k)
    sample_dir = xenium_processed_data_dir / f"{name_sample}/normalised_results/outs"

sample_normalised_counts_path = sample / f"{normalisation_method}/normalised_counts/{layer}.parquet"
sample_idx_path = sample / f"{normalisation_method}/normalised_counts/cells.parquet"

# read normalised data
X_normalised = pd.read_parquet(sample_normalised_counts_path)
X_normalised.index = pd.read_parquet(sample_idx_path).iloc[:, 0]
X_normalised.columns = X_normalised.columns.str.replace(".", "-")  # undo seurat renaming

if len(genes):
    # load raw data to reapply lower bounds QC filters
    ad = readwrite.read_xenium_sample(
        sample_dir,
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
    if segmentation == "proseg_expected":
        ad.obs_names = "proseg-" + ad.obs_names.astype(str)

    # filter cells
    ad = ad[X_normalised.index, X_normalised.columns]
    ad.layers["X_normalised"] = X_normalised
    if layer != "scale_data":  # no need to sparsify scale_data which is dense
        ad.layers["X_normalised"] = scipy.sparse.csr_matrix(ad.layers["X_normalised"])
else:
    ad = sc.AnnData(X_normalised)
    if layer != "scale_data":  # no need to sparsify scale_data which is dense
        ad.X = scipy.sparse.csr_matrix(ad.X)

print("Concatenating")
# concatenate
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
for i, lvl in enumerate(xenium_levels):
    ad.obs[lvl] = k[i]

print("Done")

# subset to genes
if len(genes):
    print("Subsetting")

    genes_found = [
        g
        for g in ad.var_names
        if (g in genes) or (g.replace(".", "-") in genes)  # possible seurat renaming
    ]

    print(f"Found {len(genes_found)} out of {len(genes)} genes.")
    ad = ad[:, genes_found].copy()
    # reapply QC to subset of genes
    preprocessing.preprocess(
        ad,
        min_counts=min_counts,
        min_genes=min_features,
        max_counts=max_counts,
        max_genes=max_features,
        min_cells=min_cells,
        save_raw=False,
    )
    # replace X
    ad.X = ad.layers["X_normalised"]


print("Computing PCA and UMAP")
# preprocess
preprocessing.preprocess(
    ad,
    normalize=False,
    log1p=False,
    scale="none",
    n_comps=n_comps,
    metric=metric,
    min_dist=min_dist,
    n_neighbors=n_neighbors,
    pca=True,
    umap=True,
    save_raw=False,
    min_counts=None,
    min_genes=None,
    max_counts=None,
    max_genes=None,
    min_cells=None,
)

# save
df_umap = pd.DataFrame(ad.obsm["X_umap"], index=ad.obs_names, columns=["UMAP1", "UMAP2"])
df_umap[xenium_levels] = ad.obs[xenium_levels]

df_umap.to_parquet(out_file)
