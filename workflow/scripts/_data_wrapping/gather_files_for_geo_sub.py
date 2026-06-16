"""
Copy XeniumRanger output files for a single sample to the GEO submission directory.

Fixed files (morphology stack, transcripts, count matrix, cell/nucleus tables) are copied
with a user-supplied prefix. 2D focus images are copied according to the major version of
XeniumRanger that produced the data, as read from the versions JSON produced by
collect_10x_versions.py:

  v1   : <root>/morphology_focus.ome.tif                         (single file)
  v2/v3: <root>/morphology_focus/morphology_focus_NNNN.ome.tif   (1 or 4 channels)
  v4+  : <root>/morphology_focus/chNNNN_<name>.ome.tif           (all channels)
"""

import re
import sys
import os
import argparse
import contextlib
import json
import shutil
from pathlib import Path


# ---------------------------------------------------------------------------
# Fixed file names relative to the XeniumRanger output root
# ---------------------------------------------------------------------------
_FIXED_FILES: list[tuple[str, str]] = [
    ("morphology.ome.tif",          "{prefix}_morphology.ome.tif"),
    ("transcripts.parquet",         "{prefix}_transcripts.parquet"),
    ("cell_feature_matrix.h5",      "{prefix}_cell_feature_matrix.h5"),
    ("cells.parquet",               "{prefix}_cells.parquet"),
    ("cell_boundaries.parquet",     "{prefix}_cell_boundaries.parquet"),
    ("nucleus_boundaries.parquet",  "{prefix}_nucleus_boundaries.parquet"),
]

# v1: single focus file at the root
_FOCUS_FILE_V1 = "morphology_focus.ome.tif"

# v2/v3 and v4: focus images live in this subdirectory
_FOCUS_DIR = "morphology_focus"

# v2/v3: morphology_focus_NNNN.ome.tif
_FOCUS_PAT_V23 = re.compile(r"^morphology_focus_\d{4}\.ome\.tif$")

# v4: chNNNN_<name>.ome.tif
_FOCUS_PAT_V4 = re.compile(r"^ch\d{4}_.+\.ome\.tif$")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tif_files_in(directory: Path) -> list[str]:
    return sorted(
        f for f in os.listdir(directory)
        if f.endswith(".ome.tif") and not f.startswith("._")
    )


def _copy(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Expected source file not found: {src}")
    shutil.copy2(src, dst)
    print(f"  {src.name}  ->  {dst.name}")


def _copy_focus_dir(focus_dir: Path, pat: re.Pattern, out_dir: Path, prefix: str) -> None:
    if not focus_dir.is_dir():
        raise FileNotFoundError(f"Focus image directory not found: {focus_dir}")
    tif_files = _tif_files_in(focus_dir)
    if not tif_files:
        raise FileNotFoundError(f"No .ome.tif files found in {focus_dir}")
    mismatched = [f for f in tif_files if not pat.match(f)]
    if mismatched:
        raise ValueError(
            f"Unexpected file names in {focus_dir} "
            f"(pattern {pat.pattern!r} not matched): {mismatched}. "
            "This may indicate a XeniumRanger major version mismatch."
        )
    for fname in tif_files:
        _copy(focus_dir / fname, out_dir / f"{prefix}_{fname}")


def gather_focus_images(root: Path, out_dir: Path, prefix: str, major_version: int) -> None:
    if major_version == 1:
        _copy(root / _FOCUS_FILE_V1, out_dir / f"{prefix}_{_FOCUS_FILE_V1}")
    elif major_version in (2, 3):
        _copy_focus_dir(root / _FOCUS_DIR, _FOCUS_PAT_V23, out_dir, prefix)
    elif major_version >= 4:
        _copy_focus_dir(root / _FOCUS_DIR, _FOCUS_PAT_V4, out_dir, prefix)
    else:
        raise ValueError(f"Unsupported XeniumRanger major version: {major_version}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy XeniumRanger files to the GEO submission directory."
    )
    parser.add_argument(
        "-i", required=True, type=str,
        help="path to the XeniumRanger output directory (raw data root)",
    )
    parser.add_argument(
        "-o", required=True, type=str,
        help="path to the output directory",
    )
    parser.add_argument(
        "--prefix", required=True, type=str,
        help="filename prefix applied to every output file",
    )
    parser.add_argument(
        "--versions", required=True, type=str,
        help="path to the versions JSON produced by collect_10x_versions.py",
    )
    parser.add_argument(
        "-l", default=None, type=str,
        help="path to the log file",
    )

    return parser.parse_args()


def _main(args: argparse.Namespace) -> None:
    if not os.path.isdir(args.i):
        raise NotADirectoryError(f"Input directory does not exist: {args.i}")
    if not os.path.isfile(args.versions):
        raise FileNotFoundError(f"Versions JSON does not exist: {args.versions}")
    os.makedirs(args.o, exist_ok=True)

    with open(args.versions, "r", encoding="utf-8") as fh:
        versions: dict = json.load(fh)

    try:
        major_version: int = int(versions["raw_data_version"]["0"])
    except (KeyError, ValueError) as exc:
        raise RuntimeError(
            f"Could not read XeniumRanger major version from {args.versions}: {exc}"
        ) from exc

    root = Path(args.i)
    out_dir = Path(args.o)

    print(f"Source          : {root}")
    print(f"Output          : {out_dir}")
    print(f"Prefix          : {args.prefix}")
    print(f"XeniumRanger v  : {major_version}.x")

    print("\nCopying fixed files...")
    for src_name, dst_template in _FIXED_FILES:
        _copy(root / src_name, out_dir / dst_template.format(prefix=args.prefix))

    print("\nCopying focus images...")
    gather_focus_images(root, out_dir, args.prefix, major_version)

    print("\nDone.")


if __name__ == "__main__":
    args = parse_args()

    log_cm = open(args.l, "w", encoding="utf-8") if args.l else contextlib.nullcontext()
    with log_cm as _log:
        if args.l:
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = _log
        try:
            _main(args)
        finally:
            if args.l:
                sys.stdout = old_stdout
                sys.stderr = old_stderr
