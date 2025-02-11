"""
This script is used to map the raw and normalised cells from the results of Proseg.
"""

import sys
import os
import argparse

import pandas as pd
from sklearn.neighbors import NearestNeighbors


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Map raw and normalised cells."
    )

    sys_args_parser.add_argument(
        "--raw",
        required=True,
        type=str,
        help="path to the cells file (cell-metadata.csv.gz) from the raw results of Proseg",
    )

    sys_args_parser.add_argument(
        "--normalised",
        required=True,
        type=str,
        help="path to the cells file (cells.parquet) from the normalised results of 10X Xenium Ranger",
    )

    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "--out",
        required=True,
        type=str,
        help="path to the output map file",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.raw):
        raise RuntimeError(
            f"Error! Input cells file from the raw results of Proseg does not exist: {ret.raw}"
        )
    if not os.path.isfile(ret.normalised):
        raise RuntimeError(
            f"Error! Input cells file from the normalised results of 10X Xenium Ranger does not exist: {ret.normalised}"
        )

    os.makedirs(os.path.dirname(os.path.abspath(ret.out)), exist_ok=True)

    return ret


def uniquify(df: pd.DataFrame) -> pd.DataFrame:
    if df.shape[0] == 1:
        return df

    ret = df.sort_values("dist")
    return ret.iloc[[0]]


def get_nearest_neighbors(
    ref_df: pd.DataFrame,
    inquery_df: pd.DataFrame,
    ref_cell_col_name: str,
    inquery_cell_col_name: str,
    *,
    x_centroid: str = "x_centroid",
    y_centroid: str = "y_centroid",
    cell_id: str = "cell_id",
) -> pd.DataFrame:
    dist, idx = (
        NearestNeighbors(n_neighbors=1)
        .fit(ref_df[[x_centroid, y_centroid]])
        .kneighbors(inquery_df[[x_centroid, y_centroid]])
    )
    ids = [x[0] for x in idx]
    dists = [x[0] for x in dist]

    return (
        pd.DataFrame(
            {
                ref_cell_col_name: ref_df[cell_id].iloc[ids].values,
                inquery_cell_col_name: inquery_df[cell_id].values,
                "dist": dists,
            }
        )
        .groupby(
            ref_cell_col_name,
        )
        .apply(
            uniquify,
        )
        .reset_index(
            drop=True,
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

    xr_cells: pd.DataFrame = pd.read_parquet(args.normalised)
    proseg_cells: pd.DataFrame = pd.read_csv(
        args.raw,
        compression="gzip",
    ).rename(
        columns={
            "cell": "cell_id",
            "centroid_x": "x_centroid",
            "centroid_y": "y_centroid",
        },
    )

    proseg2xr: pd.DataFrame = get_nearest_neighbors(
        xr_cells,
        proseg_cells,
        "xr_cell_id",
        "proseg_cell_id",
    )

    xr2proseg: pd.DataFrame = get_nearest_neighbors(
        proseg_cells,
        xr_cells,
        "proseg_cell_id",
        "xr_cell_id",
    )

    mapped = pd.merge(
        proseg2xr,
        xr2proseg,
        how="inner",
        on=["xr_cell_id", "proseg_cell_id"],
        suffixes=("_proseg2xr", "_xr2proseg"),
    )

    mapped.to_parquet(args.out, index=False)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
