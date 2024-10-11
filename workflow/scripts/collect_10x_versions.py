"""
This script is used to get the versions of 10X xeniumranger both from the raw data and from the system software.
"""

import argparse
import subprocess
import re
import json
import os
import sys
from typing import Any

from bs4 import BeautifulSoup


def parse_args():
    sys_args_parser = argparse.ArgumentParser(
        description="Collect 10X xeniumranger versions."
    )

    sys_args_parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="path to the input analysis summary file from teh raw data",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument("-o", required=True, type=str, help="output json file")

    ret = sys_args_parser.parse_args()
    if not os.path.isfile(ret.i):
        raise RuntimeError(
            f"Error! Input analysis summary file does not exist: {ret.i}"
        )

    os.makedirs(os.path.dirname(ret.o), exist_ok=True)

    return ret


def extract_raw_data_verion(file_path: str) -> str:
    with open(file_path, "r", encoding="utf-8") as fh:
        parsed = BeautifulSoup(fh, "html.parser")

    matched_target = re.search(
        r"const data = (.*)", parsed.find("script", type="text/javascript").string
    )

    if matched_target:
        try:
            software_info = json.loads(matched_target.group(1))
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Error! {e.msg}")

        if "software_table" in software_info:
            for i in software_info["software_table"]["rows"]:
                if i[0] == "Analysis version":
                    matched_ret = re.match(r"^xenium-([\d.]+)", i[1])

                    if matched_ret:
                        return matched_ret.group(1)

            raise RuntimeError(
                'Error! Cannot find related software information under "Analysis version". Please raise an issue with detailed information to the maintainer.'
            )

        raise KeyError(
            'Error! Cannot find related software information under "software_table". Please raise an issue with detailed information to the maintainer.'
        )

    raise RuntimeError(
        'Error! Cannot find related software information under "const data". Please raise an issue with detailed information to the maintainer.'
    )


def extract_system_software_version() -> str:
    matched_target = re.search(
        r"xenium-([\d.]+)",
        subprocess.check_output("xeniumranger --version", shell=True).decode("utf-8"),
    )

    if not matched_target:
        raise RuntimeError(
            "Error! Cannot find related software information from 10X xeniumranger. Please raise an issue with detailed information to the maintainer."
        )

    return matched_target.group(1)


def parse_version(version: str) -> tuple[int | str, ...]:
    ret: tuple[int | str, ...] = ()

    for i in version.split("."):
        try:
            ret = (*ret, int(i))
        except TypeError:
            ret = (*ret, i)

    return ret


def serialise_tuple(tp: tuple[Any, ...]) -> dict[str, Any]:
    ret: dict[str, Any] = {}

    for idx, val in enumerate(tp):
        ret[str(idx)] = val

    return ret


def match_versions(
    v1: tuple[int | str, ...], v2: tuple[int | str, ...]
) -> tuple[bool, ...]:
    ret: tuple[bool, ...] = ()

    for i in range(max(len(v1), len(v2))):
        if i < len(v1) and i < len(v2):
            val: bool = str(v1[i]) == str(v2[i])
        else:
            val = True

        ret = (*ret, val and (True if i == 0 else ret[-1]))

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    raw_data_version = parse_version(extract_raw_data_verion(args.i))
    system_software_version = parse_version(extract_system_software_version())

    versions: dict = {}
    versions["raw_data_version"] = serialise_tuple(raw_data_version)
    versions["system_software_version"] = serialise_tuple(system_software_version)
    versions["match"] = serialise_tuple(
        match_versions(raw_data_version, system_software_version)
    )

    with open(args.o, "w", encoding="utf-8") as fh:
        json.dump(versions, fh)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
