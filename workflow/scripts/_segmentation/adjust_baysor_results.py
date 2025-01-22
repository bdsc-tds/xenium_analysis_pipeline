"""
This script is used to adjust the results of Baysor for Xenium Ranger prior to v3.1.
"""

import argparse
import os
import sys
import re
from typing import Any
import json

import pandas as pd


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Adjust Baysor results.")

    sys_args_parser.add_argument(
        "--inseg",
        required=True,
        type=str,
        help="path to the input segmentation file",
    )
    sys_args_parser.add_argument(
        "--inpoly",
        required=True,
        type=str,
        help="path to the input polygon file",
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
        "--outpoly", required=True, type=str, help="path to the output polygon file"
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.inseg):
        raise RuntimeError(
            f"Error! Input segmentation file does not exist: {ret.inseg}"
        )
    if not os.path.isfile(ret.inpoly):
        raise RuntimeError(f"Error! Input polygon file does not exist: {ret.inpoly}")

    os.makedirs(os.path.dirname(ret.outseg), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outpoly), exist_ok=True)

    return ret


def remove_prefix2cell(polygon_info: list[dict[str, Any]]) -> set[str] | None:
    ret: set[str] = set()

    if len(polygon_info) < 1:
        return None

    for i in polygon_info:
        ret.add(i["cell"])

        matched = re.search(r".+-(\d+)$", i["cell"])

        if matched:
            i["cell"] = int(matched.group(1))
        else:
            print(f'Cannot process cell {i["cell"]}.')

    return ret


def remove_cells_in_seg(data: pd.DataFrame, cells2keep: set[str]) -> None:
    indices: set[int] = set()

    for idx, row in data.iterrows():
        if pd.isna(row.cell):
            continue

        if row.cell not in cells2keep:
            indices.add(idx)

    data.drop(index=list(indices), inplace=True)


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    with open(args.inpoly, "r", encoding="utf-8") as fh:
        polygons = json.load(fh)
    polygon_cells: set[str] | None = remove_prefix2cell(polygons["geometries"])
    assert polygon_cells is not None

    seg: pd.DataFrame = pd.read_csv(args.inseg)
    remove_cells_in_seg(seg, polygon_cells)

    with open(args.outpoly, "w", encoding="utf-8") as fh:
        json.dump(polygons, fh)
    seg.to_csv(
        args.outseg,
        index=False,
    )

    print("The results of Baysor have been processed for Xenium Ranger  prior to v3.1.")

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
