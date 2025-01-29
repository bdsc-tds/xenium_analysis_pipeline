from pathlib import Path
import argparse
import sys
import pandas as pd
import scanpy as sc

sys.path.append("workflow/scripts/utils/")
import readwrite
import preprocessing

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")
parser.add_argument("--panel", type=Path, help="Path to the panel file.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--n_neighbors", type=int, help="Number of neighbors.")
parser.add_argument("--metric", type=str, help="Distance metric to use.")
parser.add_argument("--min_dist", type=float, help="Minimum distance parameter.")
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

args = parser.parse_args()

# Access the arguments
panel = args.panel
out_file = args.out_file
n_comps = args.n_comps
n_neighbors = args.n_neighbors
metric = args.metric
min_dist = args.min_dist
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells


segmentation = panel.parents[1].stem
cohort = panel.parents[0].stem

# read xenium samples
xenium_paths = {}
for sample in (samples := panel.iterdir()):
    for replicate in (replicates := sample.iterdir()):
        k = (segmentation, cohort, panel.stem, sample.stem, replicate.stem)
        replicate_path = replicate / "normalised_results/outs"
        name = "/".join(k)

        xenium_paths[k] = replicate_path

ads = readwrite.read_xenium_samples(
    xenium_paths, anndata_only=True, transcripts=False, sample_name_as_key=False
)

# concatenate
xenium_levels = ["segmentation", "cohort", "panel", "sample", "replicate"]
for k in ads.keys():
    for i, l in enumerate(xenium_levels):
        ads[k].obs[l] = k[i]
ad_merge = sc.concat(ads)

# preprocess
preprocessing.preprocess(
    ad_merge,
    normalize=True,
    log1p=True,
    scale="none",
    n_comps=n_comps,
    metric=metric,
    min_dist=min_dist,
    n_neighbors=n_neighbors,
    pca=True,
    umap=True,
    save_raw=False,
    min_counts=min_counts,
    min_genes=min_features,
    max_counts=max_counts,
    max_genes=max_features,
    min_cells=min_cells,
)

# save
df_umap = pd.DataFrame(
    ad_merge.obsm["X_umap"], index=ad_merge.obs_names, columns=["UMAP1", "UMAP2"]
)
df_umap[xenium_levels] = ad_merge.obs[xenium_levels]

df_umap.to_parquet(out_file)
