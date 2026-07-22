"""
Convert the *new* Segger (dpeerlab/segger) output into the Baysor-style
transcript-assignment CSV + polygon GeoJSONs that `xeniumranger
import-segmentation` consumes -- exactly the contract the old
`convert_segger2baysor.py` produced, so everything downstream is unchanged.

Refer to: https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis/inputs/segmentation-inputs

Difference from the old converter
---------------------------------
The old Segger emitted `segger_transcripts.parquet` already carrying
`transcript_id`, original `cell_id` and `x_location/y_location`. The new Segger
emits `segger_segmentation.parquet` with only:
    row_index, segger_cell_id, segger_similarity, similarity_threshold, converged
where `row_index` is the 0-based positional index into the *raw* Xenium
`transcripts.parquet` (assigned BEFORE the qv>=20 / control-probe filters in
segger's preprocessor, so it maps 1:1 to the raw parquet row order). A transcript
is considered assigned only where `segger_similarity >= similarity_threshold`.

So the ONLY new work here is the input reshape: read the new parquet, keep the
above-threshold assignments, and recover `transcript_id / cell_id / x_location /
y_location` by joining `row_index` back to the raw transcripts. From that point on
the frame is identical to what the old script assumed and all downstream logic
(get_boundaries / convert_feat2geom_collection / refine_segmentation_table) is
reused verbatim.

Boundary generation: the new segger removed `segger.prediction.boundary`, so we
import `generate_boundary` from the vendored `segger_boundary.py` (same dir).
"""

import os
import sys
import re

from typing import Any

import argparse
import json
import numpy as np
import pandas as pd
import polars as pl
import geopandas as gpd

# Vendored, self-contained boundary generator (replaces the removed
# segger.prediction.boundary). Ensure this script's dir is importable regardless
# of the working directory Snakemake launches us from.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from segger_boundary import generate_boundary  # noqa: E402

CELL_ID_COLUMN_NAME = "segger_cell_id"

# New-segger output columns
COL_ROW_INDEX = "row_index"
COL_SIMILARITY = "segger_similarity"
COL_THRESHOLD = "similarity_threshold"

# Raw Xenium transcripts.parquet columns we recover (XeniumTranscriptFields)
RAW_TX_ID = "transcript_id"
RAW_CELL_ID = "cell_id"
RAW_X = "x_location"
RAW_Y = "y_location"


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Adjust *new* Segger results for Xenium Ranger import."
    )

    sys_args_parser.add_argument(
        "--inseg",
        required=True,
        type=str,
        help="path to the new segger output (segger_segmentation.parquet)",
    )
    sys_args_parser.add_argument(
        "--intx",
        required=True,
        type=str,
        help="path to the raw Xenium transcripts.parquet (to recover "
        "transcript_id / cell_id / x_location / y_location via row_index)",
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
    if not os.path.isfile(ret.intx):
        raise RuntimeError(
            f"Error! Raw transcripts file does not exist: {ret.intx}"
        )
    os.makedirs(os.path.dirname(ret.outseg), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outpolyfeat), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outpolygeom), exist_ok=True)

    return ret


def load_new_segger_assignments(inseg: str, intx: str) -> pd.DataFrame:
    """Read the new-segger parquet, keep above-threshold assignments, and recover
    transcript_id / cell_id / x_location / y_location from the raw transcripts.

    Returns a pandas DataFrame with the columns the downstream (verbatim) logic
    expects: ``segger_cell_id``, ``transcript_id``, ``cell_id``,
    ``x_location``, ``y_location``.
    """
    # --- new segger output: filter to above-threshold, non-null assignments ---
    seg = pl.read_parquet(inseg)

    # An assignment is valid only where similarity >= per-gene threshold and the
    # cell id is present (the model leaves unassigned transcripts null).
    seg = seg.filter(
        pl.col(CELL_ID_COLUMN_NAME).is_not_null()
        & pl.col(COL_SIMILARITY).is_not_null()
        & pl.col(COL_THRESHOLD).is_not_null()
        & (pl.col(COL_SIMILARITY) >= pl.col(COL_THRESHOLD))
    ).select([
        pl.col(COL_ROW_INDEX).cast(pl.Int64),
        pl.col(CELL_ID_COLUMN_NAME).cast(pl.Utf8),
    ])

    # --- raw transcripts: reproduce segger's pre-filter 0-based row index ---
    # segger's preprocessor does `scan_parquet(...).with_row_index()` BEFORE any
    # filtering, so row_index == physical raw row position. We mirror that here.
    raw = (
        pl.scan_parquet(intx, parallel="row_groups")
        .with_row_index(name=COL_ROW_INDEX)
        .select([
            pl.col(COL_ROW_INDEX).cast(pl.Int64),
            pl.col(RAW_TX_ID),
            pl.col(RAW_CELL_ID).cast(pl.Utf8),
            pl.col(RAW_X),
            pl.col(RAW_Y),
        ])
    )

    # Join from the (small) assigned side; keeps memory bounded to assignments.
    merged = (
        seg.lazy()
        .join(raw, on=COL_ROW_INDEX, how="left")
        .drop(COL_ROW_INDEX)
        .collect()
    )

    df = merged.to_pandas()
    # Downstream expects real NaN for unassigned xr cell ids. Xenium uses
    # 'UNASSIGNED' (SW >= 2.0) or '-1' (SW 1.x) as the unassigned-cell sentinel.
    df[RAW_CELL_ID] = df[RAW_CELL_ID].replace({"UNASSIGNED": np.nan, "-1": np.nan})
    return df


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
    cell_id: str = CELL_ID_COLUMN_NAME,
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
        _geo_table.set_index(CELL_ID_COLUMN_NAME),
        on=CELL_ID_COLUMN_NAME,
        how="inner",
        validate="many_to_one",
    )

    ret["is_noise"] = ret[CELL_ID_COLUMN_NAME].isna()

    ret = ret[
        ["transcript_id", "cell_id", CELL_ID_COLUMN_NAME, "dummy_cell_id", "is_noise"]
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

    # --- ONLY new step: reshape new-segger output into the old in-memory frame ---
    seg: pd.DataFrame = load_new_segger_assignments(args.inseg, args.intx)

    # --- everything below is the old convert_segger2baysor.py logic, verbatim ---
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
        f"The results of the new Segger have been converted to the format of Baysor, compatible with v0.7{'-' if args.prior2baysor07 else '+'}."
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
