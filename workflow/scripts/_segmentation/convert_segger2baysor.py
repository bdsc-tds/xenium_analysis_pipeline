"""
This script is used to adjust the results of Segger in order to import to Xenium Ranger.

Refer to: https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis/inputs/segmentation-inputs
"""

import os
import sys
import re

from typing import Any

import argparse
import json
import numpy as np
import pandas as pd
import geopandas as gpd

from segger.prediction.boundary import generate_boundary


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Adjust Segger results.")

    sys_args_parser.add_argument(
        "--inseg",
        required=True,
        type=str,
        help="path to the input segmentation file",
    )

    sys_args_parser.add_argument(
        "--prior2baysor07",
        action=argparse.BooleanOptionalAction,
        help="whether the adjusted results should be compatible with Baysor version < 0.7",
    )

    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "--outseg",
        required=True,
        type=str,
        help="path to the output segmentation file",
    )
    sys_args_parser.add_argument(
        "--outpolyfeat",
        required=True,
        type=str,
        help="path to the output polygon file with feature collection",
    )
    sys_args_parser.add_argument(
        "--outpolygeom",
        required=True,
        type=str,
        help="path to the output polygon file with geometry collection",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.inseg):
        raise RuntimeError(
            f"Error! Input segmentation file does not exist: {ret.inseg}"
        )
    os.makedirs(os.path.dirname(ret.outseg), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outpolyfeat), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outpolygeom), exist_ok=True)

    return ret


def _find_boundary_per_cell(
    df: pd.DataFrame,
    *,
    x: str = "x_location",
    y: str = "y_location",
) -> pd.DataFrame:
    try:
        geom = generate_boundary(df, x=x, y=y)
    except Exception:
        geom = None

    return pd.DataFrame(
        {
            "length": len(df),
            "geometry": [geom],
        }
    )


def get_boundaries(
    df: pd.DataFrame,
    x: str = "x_location",
    y: str = "y_location",
    cell_id: str = "segger_cell_id",
) -> gpd.GeoDataFrame:
    df = df[~df[cell_id].isna()]
    codes, levels = pd.factorize(np.unique(df[cell_id]))
    code_map = pd.DataFrame({"code": codes}, index=levels, dtype=int)
    results = (
        df.groupby(cell_id)
        .apply(_find_boundary_per_cell, include_groups=False, x=x, y=y)
        .reset_index()
        .join(code_map, on=cell_id)
    )
    results["dummy_cell_id"] = "cell-" + results["code"].astype(str)

    ret: gpd.GeoDataFrame = gpd.GeoDataFrame(
        results[[cell_id, "dummy_cell_id", "length", "geometry"]]
    )
    return ret[ret.geometry.notna() & ~ret.geometry.is_empty]


def refine_segmentation_table(
    ori_table: pd.DataFrame,
    geo_table: gpd.GeoDataFrame,
    cells2exclude: pd.DataFrame | None = None,
) -> pd.DataFrame:
    _geo_table: pd.DataFrame = pd.DataFrame(
        geo_table.drop(columns=["geometry", "length"])
    )

    if cells2exclude is not None:
        _geo_table = pd.merge(
            _geo_table,
            cells2exclude,
            on="dummy_cell_id",
            how="outer",
            indicator=True,
        )

        _geo_table = _geo_table[_geo_table["_merge"] == "left_only"].drop(
            columns=["_merge"]
        )

    ret: pd.DataFrame = ori_table.join(
        _geo_table.set_index("segger_cell_id"),
        on="segger_cell_id",
        how="inner",
        validate="many_to_one",
    )

    ret["is_noise"] = ret["segger_cell_id"].isna()

    ret = ret[
        ["transcript_id", "cell_id", "segger_cell_id", "dummy_cell_id", "is_noise"]
    ].rename(columns={"cell_id": "xr_cell_id", "dummy_cell_id": "cell"})

    return ret


def convert_feat2geom_collection(
    feat_col: dict[str, Any], prior2baysor07: bool
) -> tuple[list[str], dict[str, Any]]:
    cells2exclude: list[str] = []
    geom_col: dict[str, Any] = {"type": "GeometryCollection", "geometries": []}

    for feat in feat_col["features"]:
        if feat["geometry"]["type"] != "Polygon":
            cells2exclude.append(feat["properties"]["dummy_cell_id"])
            continue

        cell_id: str = feat["properties"]["dummy_cell_id"]
        _id: int = -1

        if prior2baysor07:
            matched = re.search(r".+-(\d+)$", cell_id)

            if matched:
                _id = int(matched.group(1))
            else:
                raise RuntimeError(f"Error! Cannot process cell {cell_id}.")

        geom_col["geometries"].append(
            {
                "type": "Polygon",
                "coordinates": feat["geometry"]["coordinates"],
                "cell": _id if prior2baysor07 and _id != -1 else cell_id,
            }
        )

    return cells2exclude, geom_col


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    seg: pd.DataFrame = pd.read_parquet(args.inseg)

    geo: gpd.GeoDataFrame = get_boundaries(seg)

    geo.to_file(args.outpolyfeat, driver="GeoJSON")

    assert os.path.isfile(args.outpolyfeat)
    with open(args.outpolyfeat, "r", encoding="utf-8") as fh:
        polygons_feat = json.load(fh)

    cells2exclude, polygons_geom = convert_feat2geom_collection(
        polygons_feat, args.prior2baysor07
    )

    refine_segmentation_table(
        seg,
        geo,
        pd.DataFrame({"dummy_cell_id": cells2exclude}),
    ).to_csv(args.outseg, index=False)

    with open(args.outpolygeom, "w", encoding="utf-8") as fh:
        json.dump(polygons_geom, fh)

    print(
        f"The results of Segger have been converted to the format of Baysor, compatible with v0.7{'-' if args.prior2baysor07 else '+'}."
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
