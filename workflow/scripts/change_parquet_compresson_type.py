"""
This script is used to change the compression type of a parquet file.
"""

import argparse
import os
import sys
from typing import Any

import pandas as pd


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Change comopression type of parquet files."
    )

    sys_args_parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="path to the input parquet file",
    )
    sys_args_parser.add_argument(
        "-t",
        type=str,
        default="snappy",
        choices=["snappy", "gzip", "brotli", "lz4", "zstd"],
        help="target compression type",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "-o", required=True, type=str, help="path to the output parquet file"
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.i):
        raise RuntimeError(f"Error! Input parquet file does not exist: {ret.i}")

    os.makedirs(os.path.dirname(ret.o), exist_ok=True)

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    table: pd.DataFrame = pd.read_parquet(args.i)
    table.to_parquet(args.o, compression=args.t)

    print(
        f'The compression type of the input parquet file has been changed to "{args.t}".'
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
