"""
This script is used to gather information from RCTD results per sample.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Gather information from RCTD results."
    )

    sys_args_parser.add_argument(
        "--in_meta",
        required=True,
        type=str,
        help="path to the input metadata file",
    )
    sys_args_parser.add_argument(
        "--in_rctd",
        required=True,
        type=str,
        help="path to the input rctd results",
    )
    sys_args_parser.add_argument(
        "--out",
        required=True,
        type=str,
        help="path to the output file",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "--normalisation_id",
        required=True,
        type=str,
        help="normalisation id",
    )
    sys_args_parser.add_argument(
        "--annotation_id",
        type=str,
        help="annotation id",
    )
    sys_args_parser.add_argument(
        "--ref_name",
        required=True,
        type=str,
        help="reference name used to run RCTD",
    )
    sys_args_parser.add_argument(
        "--ref_level",
        type=str,
        help="reference level used to run RCTD",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.in_meta):
        raise RuntimeError(f"Error! Input metadata file does not exist: {ret.in_meta}")
    if not os.path.isfile(ret.in_rctd):
        raise RuntimeError(
            f"Error! Input RCTD result file does not exist: {ret.in_rctd}"
        )

    os.makedirs(os.path.dirname(ret.out), exist_ok=True)

    return ret


def get_meta_data(
    file_path: str,
    normalisation_id: str | None,
    annotation_id: str | None,
    ref_name: str | None,
    ref_level: str | None,
) -> pd.DataFrame:
    return (
        pd.read_parquet(
            file_path,
        )
        .rename(
            columns={
                "nCount_Xenium": "ncount",
                "nFeature_Xenium": "nfeature",
            },
        )
        .filter(
            [
                "cell_id",
                "ncount",
                "nfeature",
                "segmentation_id",
                "segmentation_method",
                "sample_id",
                "condition",
                "gene_panel",
                "donor",
                "sample",
            ],
        )
        .assign(
            normalisation_id=(
                normalisation_id if normalisation_id is not None else np.nan
            ),
            annotation_id=annotation_id if annotation_id is not None else np.nan,
            reference_name=ref_name if ref_name is not None else np.nan,
            reference_level=ref_level if ref_level is not None else np.nan,
            ncell=len,
        )
    )


def get_info_from_rctd(file_path: str) -> pd.DataFrame:
    return (
        pd.read_parquet(
            file_path,
        )
        .filter(
            [
                "cell_id",
                "spot_class",
            ],
        )
        .assign(
            dummy=True,
        )
        .pivot(
            columns="spot_class",
            index="cell_id",
            values="dummy",
        )
        .fillna(
            value=False,
        )
        .reset_index(
            None,
        )
    )


def combine_meta_data_and_rctd_info(
    _meta_data: pd.DataFrame, _rctd_info: pd.DataFrame
) -> pd.DataFrame:
    return (
        pd.merge(
            _meta_data,
            _rctd_info,
            how="left",
            on="cell_id",
        )
        .fillna(value=False)
        .assign(
            unannotated=lambda x: ~(
                x["reject"]
                | x["singlet"]
                | x["doublet_certain"]
                | x["doublet_uncertain"]
            ),
            mean_ncount=lambda x: x["ncount"].mean(),
            median_ncount=lambda x: x["ncount"].median(),
            mean_nfeature=lambda x: x["nfeature"].mean(),
            median_nfeature=lambda x: x["nfeature"].median(),
        )
        .agg(
            {
                "segmentation_id": "unique",
                "segmentation_method": "unique",
                "sample_id": "unique",
                "condition": "unique",
                "gene_panel": "unique",
                "donor": "unique",
                "sample": "unique",
                "mean_ncount": "unique",
                "median_ncount": "unique",
                "mean_nfeature": "unique",
                "median_nfeature": "unique",
                "ncell": "unique",
                "reject": "sum",
                "singlet": "sum",
                "doublet_certain": "sum",
                "doublet_uncertain": "sum",
                "unannotated": "sum",
            }
        )
        .map(
            lambda x: x[0] if isinstance(x, (list, np.ndarray)) else x,
        )
        .to_frame()
        .T.assign(
            reject_prop=lambda x: x["reject"] / x["ncell"],
            singlet_prop=lambda x: x["singlet"] / x["ncell"],
            doublet_certain_prop=lambda x: x["doublet_certain"] / x["ncell"],
            doublet_uncertain_prop=lambda x: x["doublet_uncertain"] / x["ncell"],
            unannotated_prop=lambda x: (
                1.0
                - x["reject_prop"]
                - x["singlet_prop"]
                - x["doublet_certain_prop"]
                - x["doublet_uncertain_prop"]
            ),
        )
    )


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    meta_data: pd.DataFrame = get_meta_data(
        args.in_meta,
        args.normalisation_id,
        args.annotation_id,
        args.ref_name,
        args.ref_level,
    )
    rctd_info: pd.DataFrame = get_info_from_rctd(
        args.in_rctd,
    )

    combined: pd.DataFrame = combine_meta_data_and_rctd_info(
        meta_data,
        rctd_info,
    )

    combined.to_parquet(
        args.out,
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
