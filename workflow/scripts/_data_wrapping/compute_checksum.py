"""
This script is used to compute checksums of files.
"""

import sys
import os
import argparse
from pathlib import Path

from concurrent.futures import ProcessPoolExecutor

from utilnest.filesystem.validation import hash_file


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Compute file checksum.")
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
        "-t",
        type=int,
        default=1,
        help="number of threads to use (default: 1)",
    )
    sys_args_parser.add_argument(
        "--mark_algo",
        action="store_true",
    )  # default is false
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )

    me_group = sys_args_parser.add_mutually_exclusive_group(
        required=True,
    )
    me_group.add_argument(
        "--single_file",
        type=str,
        help="path to the input file for which the checksum is computed",
    )
    me_group.add_argument(
        "--batch_file",
        type=str,
        help="path to the input file containing a list of files for which the checksums are computed",
    )

    ret = sys_args_parser.parse_args()

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

    files_to_process = []
    if args.single_file:
        files_to_process.append(Path(args.single_file))
    elif args.batch_file:
        with open(args.batch_file, "r", encoding="utf-8") as f:
            files_to_process = [Path(line.strip()) for line in f if line.strip()]

    with open(args.o, "w", encoding="utf-8") as fh:
        fh.write(f"file_name,{'hash_algo,' if args.mark_algo else ''}hash_value\n")

        with ProcessPoolExecutor(max_workers=args.t) as executor:
            for _file, _hash in zip(
                files_to_process,
                executor.map(
                    hash_file,
                    files_to_process,
                    [args.algo] * len(files_to_process),
                ),
            ):
                fh.write(
                    f"{_file.name},{args.algo + ',' if args.mark_algo else ''}{_hash}\n"
                )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
