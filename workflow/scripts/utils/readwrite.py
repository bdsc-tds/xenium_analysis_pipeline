import yaml
import pandas as pd
import os
import pathlib
import json
from rds2py import read_rds
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

try:
    import msgspec
except ImportError:
    pass


def soma_to_anndata(
    soma_uri, measurement_name, X_layer_name, return_experiment=False, **kwargs
):
    """
    Export a SOMA experiment to an anndata object.

    Parameters
    ----------
    soma_uri (str): The URI (path) of the SOMA experiment.
    measurement_name (str): The measurement name (e.g., 'RNA') to extract.
    X_layer_name (str): The layer name in the X matrix to extract (e.g., 'counts').
    return_experiment (bool): Whether to return the SOMA experiment.

    Returns
    -------
    anndata.AnnData: The anndata object.
    """
    import tiledbsoma as soma
    import tiledbsoma.io

    # Open the SOMA experiment
    experiment = soma.open(soma_uri)
    # Export the SOMA experiment to anndata format
    ad = soma.io.to_anndata(
        experiment=experiment,
        measurement_name=measurement_name,
        X_layer_name=X_layer_name,
        **kwargs,
    )

    if return_experiment:
        return ad, experiment
    else:
        return ad


def xenium_specs(path):
    path = pathlib.Path(path)
    with open(path / "experiment.xenium") as f:
        specs = json.load(f)
    return specs


######### Xenium readers
def xenium_samples_files(dir_segmentation_cohort, segmentation=None, samples=None):
    """
    Get a dictionary of files for each sample in a Xenium segmentation run.

    Parameters
    ----------
    dir_segmentation_cohort (str): The directory path of the segmentation run.
    segmentation (str): The segmentation name, e.g., 'default', '10x_5um', 'baysor'.
    samples (list): The sample names to include. If None, include all samples.

    Returns
    -------
    dict: A dictionary of files for each sample.
    """
    files = {}
    for sample_path in pathlib.Path(dir_segmentation_cohort).iterdir():
        for replicate_path in sample_path.iterdir():
            sample_name = replicate_path.stem

            if samples is not None and sample_name not in samples:
                continue
            elif "corrupted" in sample_name:
                continue
            else:
                if segmentation != "default":
                    files[sample_name] = replicate_path / "normalised_results/outs"
                else:
                    files[sample_name] = replicate_path

    return files


def read_xenium_sample(
    sample_name,
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
    anndata_only=False,
):
    """
    Reads a xenium sample from a directory path.

    Parameters
    ----------
    sample_name (str): The sample name.
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
    anndata_only (bool): Whether to return only the anndata object or the full spatialdata object.

    Returns
    -------
    If anndata_only, returns a tuple of the sample name and anndata object.
    Otherwise, returns a tuple of the sample name and spatialdata object.
    """
    import spatialdata_io

    xdata = spatialdata_io.xenium(
        path,
        cells_as_circles=cells_as_circles,
        cells_boundaries=cells_boundaries,
        nucleus_boundaries=nucleus_boundaries,
        cells_labels=cells_labels,
        nucleus_labels=nucleus_labels,
        transcripts=transcripts,
        morphology_mip=morphology_mip,
        morphology_focus=morphology_focus,
        aligned_images=aligned_images,
    )

    ad = xdata["table"]
    ad.obs_names = ad.obs["cell_id"].values

    metrics_path = pathlib.Path(path) / "metrics_summary.csv"
    if metrics_path.exists():
        ad.uns["metrics_summary"] = pd.read_csv(metrics_path)
    else:
        print("metrics_summary.csv not found at:", metrics_path)

    if anndata_only:
        return sample_name, ad
    return sample_name, xdata


def read_xenium_samples(
    data_dirs,
    cells_as_circles=False,
    cells_boundaries=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    anndata_only=False,
    sample_name_as_key=True,
):
    """
    Reads in a dictionary of sample directories and returns a dictionary of
    AnnData objects or spatialdata objects depending on the anndata_only flag.

    Parameters
    ----------
    data_dirs : dict or list
        A dictionary of sample directories or a list of paths to sample directories.
    cells_as_circles : bool, optional
        Whether to include cell boundary data as circles, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary data, by default False
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
    anndata_only : bool, optional
        Whether to only return an AnnData object, by default False
    sample_name_as_key: bool, optional
        Whether to use the sample name as the key in the return dictionary, otherwise returns full path

    Returns
    -------
    dict
        A dictionary of sample names mapped to AnnData objects or spatialdata objects.
    """
    if isinstance(data_dirs, list):
        sample_names = [
            pathlib.Path(path).stem if sample_name_as_key else path
            for path in data_dirs
        ]
        data_dirs = {
            sample_name: path for sample_name, path in zip(sample_names, data_dirs)
        }

    # Parallel processing
    xdatas = {}
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                read_xenium_sample,
                sample_name,
                path,
                cells_as_circles,
                cells_boundaries,
                nucleus_boundaries,
                cells_labels,
                nucleus_labels,
                transcripts,
                morphology_mip,
                morphology_focus,
                aligned_images,
                anndata_only,
            )
            for sample_name, path in data_dirs.items()
        ]

        for future in as_completed(futures):
            try:
                sample_name, result = future.result()
                xdatas[sample_name] = result
            except Exception as e:
                print(f"Error processing {e}")

    return xdatas


######### RCTD readers


def read_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def read_json_msgspec(file_path):
    with open(file_path, "rb") as file:
        return msgspec.json.decode(file.read())


def _rds2py_dict_to_df(r_obj_df, mode="results_df"):
    if mode == "results_df":
        r_obj_df_columns = r_obj_df["attributes"]["names"]["data"]
        r_obj_df_index = r_obj_df["attributes"]["row.names"]["data"]
        pandas_df = pd.DataFrame(
            [r_obj_df["data"][i]["data"] for i in range(len(r_obj_df["data"]))],
            index=r_obj_df_columns,
            columns=r_obj_df_index,
        ).T
    elif mode == "weights":
        r_obj_df_columns = r_obj_df["attributes"]["dimnames"]["data"][1]["data"]
        r_obj_df_index = r_obj_df["attributes"]["dimnames"]["data"][0]["data"]
        r_obj_df["data"] = r_obj_df["data"].reshape(
            r_obj_df["attributes"]["dim"]["data"], order="F"
        )
        pandas_df = pd.DataFrame(
            r_obj_df["data"], index=r_obj_df_index, columns=r_obj_df_columns
        )

    return pandas_df


def read_rctd_sample(sample_name, rctd_results_path):
    """
    Reads RCTD results from a single sample and returns a dictionary containing:

    - results_df: a pandas DataFrame with columns to be added to the anndata object's obs
    - weights: a pandas Series with the weights for each cell for the given reference
    - weights_doublet: (not implemented) a pandas Series with the weights for each cell for doublets for the given reference
    - singlet_scores: (not implemented) a pandas Series with the singlet scores for each cell for the given reference

    Parameters
    ----------
    sample_name: str
        The name of the sample
    rctd_results_path : str
        The path to the sample RCTD results
    rsuffix : str, optional
        The suffix to append to the reference name when storing to the anndata objects

    Returns
    -------
    A tuple containing the sample name and the results dictionary
    """

    r_obj = read_rds(rctd_results_path)

    results = r_obj["attributes"]["results"]
    results_keys = results["attributes"]["names"]["data"]
    results_keys_idx = {k: results_keys.index(k) for k in results_keys}

    pandas_results = {}
    for k in ["results_df", "weights"]:
        pandas_results[k] = _rds2py_dict_to_df(
            results["data"][results_keys_idx[k]], mode=k
        )

    return sample_name, pandas_results


def read_rctd_samples(ads, rctd_results_paths, prefix=""):
    """
    Read RCTD results into anndata objects in parallel using ProcessPoolExecutor.

    Parameters
    ----------
    ads : dict of anndata.AnnData
        The anndata objects to be updated.
    rctd_results_paths : str
        The directory containing the RCTD results.
    prefix : str, optional
        The prefix to append to the reference name when storing to the anndata objects.

    Returns
    -------
    None
    """

    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(
                read_rctd_sample,
                sample_name,
                rctd_results_paths[sample_name],
            ): sample_name
            for sample_name in ads.keys()
        }

        # Update anndata objects in the parent process
        for future in as_completed(futures):
            try:
                sample_name, results = future.result()
                if results:
                    ad = ads[sample_name]
                    ad.obs = ad.obs.join(results["results_df"].add_prefix(prefix))
                    ad.uns[f"{prefix}_weights"] = results["weights"]

            except Exception as e:
                print(f"Error processing sample {futures[future]}: {e}")


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
