import scanpy as sc
import torch
import os
import pathlib
import pandas as pd
import numpy as np
import importlib
import sklearn
import scipy
import seaborn as sns
import warnings


if importlib.util.find_spec("rapids_singlecell") is not None:
    # Import gpu libraries, Initialize rmm and cupy
    import rapids_singlecell as rsc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.disable_logging()


def check_gpu_availability():
    import torch

    if torch.cuda.is_available():
        return True
    else:
        return False


def split_batches(adata, batch, hvg=None, return_categories=False):
    """Split batches and preserve category information (taken from https://github.com/theislab/scib)

    :param adata:
    :param batch: name of column in ``adata.obs``. The data type of the column must be of ``Category``.
    :param hvg: list of highly variable genes
    :param return_categories: whether to return the categories object of ``batch``
    """
    split = []
    batch_categories = adata.obs[batch].cat.categories
    if hvg is not None:
        adata = adata[:, hvg]
    for i in batch_categories:
        split.append(adata[adata.obs[batch] == i].copy())
    if return_categories:
        return split, batch_categories
    return split


def merge_adata(*adata_list, **kwargs):
    """Merge adatas from list while removing duplicated ``obs`` and ``var`` columns

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """

    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(
                    f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`."
                )
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return sc.AnnData.concatenate(*adata_list, **kwargs)


def scale_batch(adata, batch):
    """Batch-aware scaling of count matrix (taken from https://github.com/theislab/scib)

    Scaling counts to a mean of 0 and standard deviation of 1 using ``scanpy.pp.scale`` for each batch separately.

    :param adata: ``anndata`` object with normalised and log-transformed counts
    :param batch: ``adata.obs`` column
    """

    # Store layers for after merge (avoids vstack error in merge)
    adata_copy = adata.copy()
    tmp = dict()
    for lay in list(adata_copy.layers):
        tmp[lay] = adata_copy.layers[lay]
        del adata_copy.layers[lay]

    split = split_batches(adata_copy, batch)

    for i in split:
        sc.pp.scale(i)

    adata_scaled = merge_adata(*split, batch_key=batch, index_unique=None)
    # Reorder to original obs_name ordering
    adata_scaled = adata_scaled[adata.obs_names]

    # Add layers again
    for key in tmp:
        adata_scaled.layers[key] = tmp[key]

    del tmp
    del adata_copy

    return adata_scaled


def preprocess(
    adata,
    batch_key="dataset_merge_id",
    normalize=False,  # Normalize total counts
    log1p=False,  # Log1p transform
    pca=False,  # Perform PCA
    scale="none",  # Scale data
    umap=False,  # Perform UMAP
    n_comps=50,  # Number of PCA components
    n_neighbors=30,  # Number of neighbors for kNN
    metric="cosine",  # Metric for kNN
    min_dist=0.3,
    backend="gpu",  # "gpu" or "cpu"
    device=0,  # Device ID for GPU backend
    save_raw=True,  # Whether to save raw data
    verbose=True,  # Optional verbose output
    min_counts=None,
    min_genes=None,
    max_counts=None,
    max_genes=None,
    min_cells=None,
):
    """
    Preprocess anndata object.

    Parameters
    ----------
    adata
        anndata object
    batch_key
        column name in adata.obs to use for batch information
    normalize
        whether to normalize total counts
    log1p
        whether to apply log1p transformation
    pca
        whether to perform PCA
    scale
        whether to scale data, can be "all", "batch" or "none"
    umap
        whether to perform UMAP
    n_comps
        number of PCA components
    n_neighbors
        number of neighbors for kNN
    metric
        metric for kNN
    min_dist
        minimum distance for UMAP
    backend
        whether to use "gpu" or "cpu" backend
    device
        device ID for GPU backend
    save_raw
        whether to save raw data in layers['counts']
    verbose
        whether to print verbose output
    min_counts
        minimum number of counts to filter cells
    min_genes
        minimum number of genes to filter cells
    max_counts
        maximum number of counts to filter cells
    max_genes
        maximum number of genes to filter cells
    min_cells
        minimum number of cells to filter genes

    Returns
    -------
    adata
        preprocessed anndata object
    """
    if "preprocess" in adata.uns:
        print("Warning: preprocess key already found in adata.uns")

    if save_raw:
        if verbose:
            print("Saving raw counts in layers['counts']...")
        adata.layers["counts"] = adata.X

    n_cells_raw = adata.shape[0]
    n_genes_raw = adata.shape[1]

    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if min_counts is not None:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    if max_genes is not None:
        sc.pp.filter_cells(adata, max_genes=max_genes)

    if verbose:
        print("Removed", n_cells_raw - adata.shape[0], " cells...")
        print("Removed", n_genes_raw - adata.shape[1], " genes...")

    ### optional switch to GPU backend ###
    if backend == "gpu":
        if not check_gpu_availability():
            print("GPU not available. Switching to CPU backend...")
            xsc = sc
            backend = "cpu"
        else:
            # allow memory oversubscription, transfer data to GPU
            rmm.reinitialize(managed_memory=True, devices=device)
            cp.cuda.set_allocator(rmm_cupy_allocator)
            if verbose:
                print("Transferring data to GPU...")
            rsc.get.anndata_to_GPU(adata)
            xsc = rsc
    else:
        xsc = sc

    ### preprocessing ###
    if normalize:
        if verbose:
            print("Normalizing total counts...")
        xsc.pp.normalize_total(adata, target_sum=1e4)
    if log1p:
        if verbose:
            print("Applying log1p transformation...")
        xsc.pp.log1p(adata)
    if scale == "all":
        if verbose:
            print("Scaling data...")
        adata = scale_batch(adata, batch_key)
    elif scale == "batch":
        if verbose:
            print("Batch-aware scaling of data...")
        xsc.pp.scale(adata)
    if pca:
        if verbose:
            print("Performing PCA...")
        xsc.tl.pca(adata, n_comps=n_comps)
    if umap:
        if verbose:
            print("Performing UMAP...")
        xsc.pp.neighbors(adata, n_pcs=n_comps, n_neighbors=n_neighbors, metric=metric)
        xsc.tl.umap(adata, min_dist=min_dist)

    # Transfer data back to CPU if using GPU backend
    if backend == "gpu":
        if verbose:
            print("Transferring data back to CPU...")
        rsc.get.anndata_to_CPU(adata)
    adata.uns["preprocess"] = dict(
        normalize=normalize,
        log1p=log1p,
        pca=pca,
        scale=scale,
        umap=umap,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        metric=metric,
        backend=backend,
        device=device,
    )
