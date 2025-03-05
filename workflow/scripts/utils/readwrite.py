import dask

dask.config.set({"dataframe.query-planning": False})

import pandas as pd
import json
import h5py
import numpy as np
import scipy
import geopandas as gpd
import dask.dataframe as dd
import spatialdata
import spatialdata_io
import warnings
import anndata as ad
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from types import MappingProxyType
from spatialdata.models import (
    ShapesModel,
    PointsModel,
)


######### Xenium readers


def read_xenium_specs(xenium_specs_file):
    xenium_specs_file = Path(xenium_specs_file)
    with open(xenium_specs_file) as f:
        specs = json.load(f)
    return specs


def xenium_proseg(
    path,
    cells_boundaries=True,
    cells_boundaries_layers=True,
    nucleus_boundaries=False,
    cells_as_circles=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=True,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    n_jobs=1,
    imread_kwargs=MappingProxyType({}),
    image_models_kwargs=MappingProxyType({}),
    labels_models_kwargs=MappingProxyType({}),
    cells_metadata=True,
    xeniumranger_dir=None,
    xenium_specs=True,
    pandas_engine="pyarrow",
    verbose=False,
):
    """
    Reads a Xenium segmentation run and returns a `SpatialData` object.

    Parameters
    ----------
    path : str
        The directory path of the segmentation run.
    cells_boundaries : bool or str, optional
        The file path of the cell boundaries GeoJSON file.
    cells_boundaries_layers : bool or str, optional
        The file path of the cell boundaries layers GeoJSON file.
    nucleus_boundaries : bool or str, optional
        Whether to read nucleus boundaries. Not implemented.
    cells_as_circles : bool or str, optional
        Whether to read cells as circles.
    cells_labels : bool or str, optional
        Whether to read cells labels. Not implemented.
    nucleus_labels : bool or str, optional
        Whether to read nucleus labels. Not implemented.
    transcripts : bool or str, optional
        The file path of the transcripts CSV file.
    morphology_mip : bool, optional
        Whether to read morphology MIP images.
    morphology_focus : bool, optional
        Whether to read morphology focus images.
    aligned_images : bool, optional
        Whether to read aligned images.
    cells_table : bool or str, optional
        The file path of the cells table CSV file.
    n_jobs : int, optional
        The number of jobs to use. Not implemented.
    imread_kwargs : dict, optional
        Keyword arguments to pass to `imread`.
    image_models_kwargs : dict, optional
        Keyword arguments to pass to `ImageModel`.
    labels_models_kwargs : dict, optional
        Keyword arguments to pass to `LabelsModel`.
    cells_metadata : bool or str, optional
        The file path of the cells metadata CSV file.
    xeniumranger_dir : str, optional
        The directory path of the XeniumRanger run.
    xenium_specs : bool or str, optional
        The file path of the Xenium specs file.
    pandas_engine : str, optional
        The pandas engine to use when reading CSV files.

    Returns
    -------
    sdata : SpatialData
        The `SpatialData` object.
    """
    path = Path(path)

    # unsupported options compared to spatialdata_io.xenium
    if nucleus_boundaries:
        raise ValueError("reading nucleus_boundaries not implemented for proseg")
    if cells_labels:
        raise ValueError("reading cells_labels not implemented for proseg")
    if nucleus_labels:
        raise ValueError("reading nucleus_labels not implemented for proseg")
    if n_jobs > 1:
        raise ValueError("n_jobs>1 not supported")

    # default expected file paths
    def parse_arg(arg, default):
        if arg:
            return default
        elif isinstance(arg, str):
            return Path(arg)
        else:
            return arg

    cells_metadata = parse_arg(cells_metadata, path / "cell-metadata.csv.gz")
    cells_boundaries = parse_arg(cells_boundaries, path / "cell-polygons.geojson.gz")
    cells_boundaries_layers = parse_arg(
        cells_boundaries_layers, path / "cell-polygons-layers.geojson.gz"
    )
    transcripts = parse_arg(transcripts, path / "transcript-metadata.csv.gz")
    cells_table = parse_arg(cells_table, path / "expected-counts.csv.gz")
    xeniumranger_dir = parse_arg(xeniumranger_dir, None)
    if xeniumranger_dir is not None or isinstance(xenium_specs, str):
        xenium_specs = parse_arg(xenium_specs, xeniumranger_dir / "experiment.xenium")

    ### images
    if morphology_mip or morphology_focus or aligned_images:
        if verbose:
            print("Reading images...")
        sdata_images = spatialdata_io.xenium(
            xeniumranger_dir,
            cells_table=False,
            cells_as_circles=False,
            cells_boundaries=False,
            nucleus_boundaries=False,
            cells_labels=False,
            nucleus_labels=False,
            transcripts=False,
            morphology_mip=morphology_mip,
            morphology_focus=morphology_focus,
            aligned_images=aligned_images,
            imread_kwargs=imread_kwargs,
            image_models_kwargs=image_models_kwargs,
            labels_models_kwargs=labels_models_kwargs,
        )

        images = sdata_images.images
    else:
        images = {}

    ### tables
    region = "cell_polygons"
    region_key = "region"
    instance_key = "cell_id"

    # flag columns not corresponding to genes
    if isinstance(cells_table, Path):
        if verbose:
            print("Reading cells table...")
        df_table = pd.read_csv(cells_table, engine=pandas_engine)

        control_columns = df_table.columns.str.contains(
            "|".join(["BLANK_", "Codeword", "NegControl"])
        )

        table = ad.AnnData(
            df_table.iloc[:, ~control_columns],
            uns={
                "spatialdata_attrs": {
                    "region": region,
                    "region_key": region_key,
                    "instance_key": instance_key,
                }
            },
        )

        if isinstance(cells_metadata, Path):
            if verbose:
                print("Reading cells metadata...")

            df_cells_metadata = pd.read_csv(
                cells_metadata, engine=pandas_engine
            ).rename(columns={"cell": "cell_id"})
            table.obs = pd.concat(
                (df_cells_metadata, df_table.iloc[:, control_columns]), axis=1
            )
            table.obsm["spatial"] = table.obs[["centroid_x", "centroid_y"]].values
            table.obs = table.obs
        table.obs[region_key] = region

        # sparsify .X
        table.X = scipy.sparse.csr_matrix(table.X)
        tables = {"table": table}
    else:
        tables = {}

    ### labels
    # not implemented
    labels = {}

    ### points
    if isinstance(transcripts, Path):
        if verbose:
            print("Reading transcripts...")

        df_transcripts = dd.read_csv(transcripts, blocksize=None).rename(
            columns={"gene": "feature_name", "assignment": "cell_id"}
        )
        points = {"transcripts": PointsModel.parse(df_transcripts)}
    else:
        points = {}

    ### shapes
    shapes = {}

    # read specs
    if isinstance(xenium_specs, Path):
        if verbose:
            print("Reading specs...")
        specs = read_xenium_specs(xenium_specs)

        # get xenium pixel size
        scale = spatialdata.transformations.Scale(
            [1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y")
        )
        transformations = {"global": scale}
    else:
        transformations = None
        if isinstance(cells_boundaries, Path) or isinstance(
            cells_boundaries_layers, Path
        ):
            warnings.warn(
                """
                Couldn't load xenium specs file with pixel size. 
                Not applying scale transformations to shapes.
                Please specify xeniumranger_dir or xenium_specs
                """
            )

    # read cells boundaries
    if isinstance(cells_boundaries, Path):
        if verbose:
            print("Reading cells boundaries...")

        df_cells_boundaries = gpd.read_file(
            "gzip://" + cells_boundaries.as_posix()
        ).rename(columns={"cell": "cell_id"})
        shapes["cells_boundaries"] = ShapesModel.parse(
            df_cells_boundaries, transformations=transformations
        )

    # read cells boundaries layers
    if isinstance(cells_boundaries_layers, Path):
        if verbose:
            print("Reading cells boundaries layers...")

        df_cells_boundaries_layers = gpd.read_file(
            "gzip://" + cells_boundaries_layers.as_posix()
        ).rename(columns={"cell": "cell_id"})
        shapes["cells_boundaries_layers"] = ShapesModel.parse(
            df_cells_boundaries_layers, transformations=transformations
        )

    # convert cells boundaries to circles
    if cells_as_circles:
        if verbose:
            print("Converting cells boundaries to circle...")

        shapes["cells_boundaries_circles"] = spatialdata.to_circles(
            shapes["cells_boundaries"]
        )

    ### sdata
    sdata = spatialdata.SpatialData(
        images=images, labels=labels, points=points, shapes=shapes, tables=tables
    )

    return sdata


def read_xenium_sample(
    path,
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
    cells_table=True,
    anndata=False,
    xeniumranger_dir=None,
    sample_name=None,
):
    """
    Reads a xenium sample from a directory path.

    Parameters
    ----------
    path (str): The directory path of the segmentation run.
    cells_as_circles (bool): Whether to include cell polygons as circles or not.
    cells_boundaries (bool): Whether to include cell boundaries or not.
    nucleus_boundaries (bool): Whether to include nucleus boundaries or not.
    cells_labels (bool): Whether to include cell labels or not.
    nucleus_labels (bool): Whether to include nucleus labels or not.
    transcripts (bool): Whether to include transcript locations or not.
    morphology_mip (bool): Whether to include morphology MIP or not.
    morphology_focus (bool): Whether to include morphology focus or not.
    aligned_images (bool): Whether to include aligned images or not.
    cells_table (bool): Whether to include cells table or not.
    anndata (bool): Whether to return only the anndata object or the full spatialdata object.
    xeniumranger_dir (str): Path to xeniumranger output directory (for proseg raw only)
    sample_name (str): The sample name.

    Returns
    -------
    If anndata, returns a tuple of the sample name and anndata object.
    Otherwise, returns a tuple of the sample name and spatialdata object.

    If sample_name is None, sample_name is not returned
    """
    path = Path(path)
    kwargs = dict(
        cells_as_circles=cells_as_circles,
        cells_boundaries=cells_boundaries,
        nucleus_boundaries=nucleus_boundaries,
        cells_labels=cells_labels,
        nucleus_labels=nucleus_labels,
        transcripts=transcripts,
        morphology_mip=morphology_mip,
        morphology_focus=morphology_focus,
        aligned_images=aligned_images,
        cells_table=cells_table,
    )

    # automatically check whether path is a folder with proseg raw outputs or in xeniumranger format
    if (path / "expected-counts.csv.gz").exists():
        reader = xenium_proseg
        kwargs["cells_boundaries_layers"] = cells_boundaries_layers
        kwargs["xeniumranger_dir"] = xeniumranger_dir
    else:
        reader = spatialdata_io.xenium

    sdata = reader(path, **kwargs)

    adata = sdata["table"]
    adata.obs_names = adata.obs["cell_id"].values

    metrics_path = path / "metrics_summary.csv"
    if metrics_path.exists():
        adata.uns["metrics_summary"] = pd.read_csv(metrics_path)
    else:
        print("metrics_summary.csv not found at:", metrics_path)

    if sample_name is None:
        if anndata:
            return adata
        else:
            return sdata
    else:
        if anndata:
            return sample_name, adata
        else:
            return sample_name, sdata


def read_xenium_samples(
    data_dirs,
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
    cells_table=True,
    anndata=False,
    sample_name_as_key=True,
    xeniumranger_dir=None,
):
    """
    Reads in a dictionary of sample directories and returns a dictionary of
    AnnData objects or spatialdata objects depending on the anndata flag.

    Parameters
    ----------
    data_dirs : dict or list
        A dictionary of sample directories or a list of paths to sample directories.
    cells_as_circles : bool, optional
        Whether to include cell boundary data as circles, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary data, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary layers data (for proseg raw only), by default False
    nucleus_boundaries : bool, optional
        Whether to include nucleus boundary data, by default False
    cells_labels : bool, optional
        Whether to include cell labels, by default False
    nucleus_labels : bool, optional
        Whether to include nucleus labels, by default False
    transcripts : bool, optional
        Whether to include transcript data, by default False
    morphology_mip : bool, optional
        Whether to include morphology data at the maximum intensity projection, by default False
    morphology_focus : bool, optional
        Whether to include morphology data at the focus, by default False
    aligned_images : bool, optional
        Whether to include aligned images, by default False
    cells_table (bool):
        Whether to include cells table or not, by default True
    anndata : bool, optional
        Whether to only return an AnnData object, by default False
    sample_name_as_key: bool, optional
        Whether to use the sample name as the key in the return dictionary,
        otherwise returns full path as key
    xeniumranger_dir: str, optional
        Path to xeniumranger output dir (for proseg raw only)

    Returns
    -------
    dict
        A dictionary of sample names mapped to AnnData objects or spatialdata objects.
    """
    if isinstance(data_dirs, list):
        sample_names = [
            Path(path).stem if sample_name_as_key else path for path in data_dirs
        ]
        data_dirs = {
            sample_name: path for sample_name, path in zip(sample_names, data_dirs)
        }

    # Parallel processing
    sdatas = {}
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                read_xenium_sample,
                path,
                cells_as_circles,
                cells_boundaries,
                cells_boundaries_layers,
                nucleus_boundaries,
                cells_labels,
                nucleus_labels,
                transcripts,
                morphology_mip,
                morphology_focus,
                aligned_images,
                anndata,
                xeniumranger_dir,
                sample_name,
            )
            for sample_name, path in data_dirs.items()
        ]

        for future in as_completed(futures):
            try:
                sample_name, result = future.result()
                sdatas[sample_name] = result
            except Exception as e:
                print(f"Error processing {e}")

    return sdatas


###### coexpression files readers
def read_coexpression_file(k, method, target_count, results_dir):
    """
    Worker function to read the coexpression and positivity rate parquet for a single sample.

    Parameters
    ----------
    k : tuple
        The sample name tuple (segmentation, cohort, panel, sample, replicate).
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    results_dir : Path
        The directory containing the coexpression results.

    Returns
    -------
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    cc : pd.DataFrame
        The coexpression matrix.
    pos_rate : pd.Series
        The positivity rate.
    """
    out_file_coexpr = (
        results_dir
        / f"coexpression/{'/'.join(k)}/coexpression_{method}_{target_count}.parquet"
    )
    out_file_pos_rate = (
        results_dir
        / f"coexpression/{'/'.join(k)}/positivity_rate_{method}_{target_count}.parquet"
    )

    cc = pd.read_parquet(out_file_coexpr)
    pos_rate = pd.read_parquet(out_file_pos_rate)[0]
    return method, target_count, cc, pos_rate


def read_coexpression_files(cc_paths, results_dir):
    """
    Reads coexpression parquet files for multiple methods and target counts in parallel using ThreadPoolExecutor.

    Parameters
    ----------
    cc_paths : list of tuples
        A list of tuples containing the key `k`, i.e., a sample name tuple (segmentation, cohort, panel, sample, replicate)
        and the method and target count to read.
    results_dir : str
        The directory containing the coexpression results.

    Returns
    -------
    CC : dict
        A dictionary with the coexpression matrices for each method and target count.
    pos_rate : dict
        A dictionary with the positivity rates for each method and target count.
    """
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                read_coexpression_file, k, method, target_count, results_dir
            )
            for k, method, target_count in cc_paths
        ]

        CC = {}
        pos_rate = {}
        for future in as_completed(futures):
            method, target_count, cc, pr = future.result()
            k = cc_paths[futures.index(future)][
                0
            ]  # Retrieve the `k` corresponding to this future

            if k not in CC:
                CC[k] = {}
            if k not in pos_rate:
                pos_rate[k] = {}

            CC[k][method, target_count] = cc
            pos_rate[k][method, target_count] = pr
    return CC, pos_rate


######### 10x writers


def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    taken from https://github.com/scverse/anndata/issues/595

    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    if ".h5" not in file:
        file = f"{file}.h5"
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if "genome" not in adata.var:
        adata.var["genome"] = "undefined"
    if "feature_types" not in adata.var:
        adata.var["feature_types"] = "Gene Expression"
    if "gene_ids" not in adata.var:
        adata.var["gene_ids"] = adata.var_names

    w = h5py.File(file, "w")
    grp = w.create_group("matrix")
    grp.create_dataset(
        "barcodes",
        data=np.array(adata.obs_names, dtype=f"|S{str_max(adata.obs_names)}"),
    )
    grp.create_dataset(
        "data", data=np.array(adata.X.data, dtype=f"<i{int_max(adata.X.data)}")
    )
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset(
        "feature_type",
        data=np.array(
            adata.var.feature_types, dtype=f"|S{str_max(adata.var.feature_types)}"
        ),
    )
    ftrs.create_dataset(
        "genome",
        data=np.array(adata.var.genome, dtype=f"|S{str_max(adata.var.genome)}"),
    )
    ftrs.create_dataset(
        "id",
        data=np.array(adata.var.gene_ids, dtype=f"|S{str_max(adata.var.gene_ids)}"),
    )
    ftrs.create_dataset(
        "name", data=np.array(adata.var.index, dtype=f"|S{str_max(adata.var.index)}")
    )
    grp.create_dataset(
        "indices", data=np.array(adata.X.indices, dtype=f"<i{int_max(adata.X.indices)}")
    )
    grp.create_dataset(
        "indptr", data=np.array(adata.X.indptr, dtype=f"<i{int_max(adata.X.indptr)}")
    )
    grp.create_dataset(
        "shape",
        data=np.array(list(adata.X.shape)[::-1], dtype=f"<i{int_max(adata.X.shape)}"),
    )


######### Add metadata to anndata
def add_metadata2ad(
    adata: ad.AnnData,
    sample_id: str,
    segmentation_id: str,
    condition: str,
    gene_panel: str,
    donor: str,
    sample: str,
    segmentation_method: str,
) -> ad.AnnData:
    """
    Adds metadata to an AnnData object.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to add metadata to.
    sample_id : str
        The sample ID.
    segmentation_id : str
        The segmentation ID.
    condition : str
        The condition.
    gene_panel : str
        The gene panel.
    donor : str
        The donor.
    sample : str
        The sample.
    segmentation_method : str
        The segmentation method.

    Returns
    -------
    ad.AnnData
        The AnnData object with metadata added.
    """
    adata.obs["sample_id"] = sample_id
    adata.obs["segmentation_id"] = segmentation_id
    adata.obs["condition"] = condition
    adata.obs["gene_panel"] = gene_panel
    adata.obs["donor"] = donor
    adata.obs["sample"] = sample
    adata.obs["segmentation_method"] = segmentation_method

    return adata
