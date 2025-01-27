"""
This script is used to clean the output directory of Segger's prediction step.
"""

import argparse
import os
import sys
import re
import shutil

import pyarrow as pa
import pyarrow.parquet as pq

COMPRESSED_PARQUET_PREFIX: str = "compressed_"

SUBDIR_PAT: str = r"^segger_embedding_.+"

output_file_names: list[str] = [
    "transcripts_df.parquet",
    "segger_adata.h5ad",
    "segger_transcripts.parquet",
    "segmentation_log.json",
]


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Clean the output directory of Segger's prediction step."
    )

    sys_args_parser.add_argument(
        "--dir",
        required=True,
        type=str,
        help="path to the output directory of Segger's prediction step",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isdir(ret.dir):
        raise RuntimeError(
            f"Error! The output directory of Segger's prediction step does not exist: {ret.dir}"
        )

    return ret


def process_outputs(path: str):
    subdir: str = ""
    existing_file_names: list[str] = []

    for dir_path, dir_names, file_names in os.walk(path):
        subdirs: list[str] = [
            i for i in dir_names if re.match(SUBDIR_PAT, i, flags=re.IGNORECASE)
        ]

        for i in file_names:
            if i in output_file_names:
                existing_file_names.append(i)

        if len(subdirs) > 1:
            print(
                f"Warning! Expecting one subdirectory in {path} but finding {len(subdirs)}: {",".join(subdirs)}. Taking the most recently modified one."
            )
            for i in subdirs:
                if subdir == "":
                    subdir = i
                else:
                    if os.path.getmtime(os.path.join(dir_path, i)) > os.path.getmtime(
                        os.path.join(dir_path, subdir)
                    ):
                        subdir = i
                    elif os.path.getmtime(
                        os.path.join(dir_path, i)
                    ) == os.path.getmtime(os.path.join(dir_path, subdir)):
                        raise RuntimeError(
                            f"Error! Two subdirectories have the same modified time: {subdir} and {i}."
                        )
        elif len(subdirs) == 1:
            subdir = subdirs[0]
        else:
            if len(existing_file_names) == len(output_file_names):
                print(f"All output files have been found in {path}. Nothing to do.")
                return None

        subdir = os.path.join(dir_path, subdir)
        break

    for dir_path, dir_names, file_names in os.walk(subdir):
        if not all(i in file_names + dir_names for i in output_file_names):
            raise FileNotFoundError(
                f"Error! Not all output files have been found in {subdir}. Missing files: {','.join([i for i in output_file_names if i not in file_names and i not in dir_names])}"
            )

        for i in file_names:
            shutil.move(os.path.join(dir_path, i), path)

        for i in dir_names:
            shutil.move(compress_parquet_files(dir_path, i), path)

        break

    shutil.rmtree(subdir, ignore_errors=True)
    return None


def compress_parquet_files(dir_path: str, dir_name: str) -> str:
    _dir: str = os.path.join(dir_path, dir_name)
    assert os.path.exists(_dir)

    tables = [
        pq.read_table(os.path.join(_dir, file))
        for file in os.listdir(_dir)
        if file.endswith(".parquet")
    ]
    table = pa.concat_tables(tables)

    output_file_path: str = os.path.join(
        dir_path,
        (
            COMPRESSED_PARQUET_PREFIX + dir_name + ""
            if dir_name.endswith(".parquet")
            else ".parquet"
        ),
    )
    pq.write_table(table, output_file_path)

    shutil.rmtree(_dir, ignore_errors=True)
    return output_file_path


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    process_outputs(args.dir)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
