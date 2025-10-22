"""
This script is used to compute checksums of files.
"""

import sys
import os
import argparse
from pathlib import Path

from utilnest.filesystem.validation import hash_file


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Compute file checksum.")
    sys_args_parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="path to the input file",
    )
    sys_args_parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="path to the output checksum file",
    )
    sys_args_parser.add_argument(
        "--algo",
        type=str,
        default="sha512",
        choices=["md5", "sha256", "sha512"],
        help="checksum algorithm to use (default: sha512)",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.i):
        raise RuntimeError(f"Error! Input file does not exist: {ret.i}")

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

    checksum = hash_file(
        args.i,
        algo=args.algo,
    )

    with open(args.o, "w", encoding="utf-8") as f:
        f.write("file_name,hash_algo,checksum\n")
        f.write(f"{Path(args.i).name},{args.algo},{checksum}\n")

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
