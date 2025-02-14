"""
Compute coexpression with given method and target counts per segmentation method per sample.
"""

import os
import argparse
import sys

import spatialdata_io

import coexpression as ce

sys.path.append("workflow/scripts/utils/readwrite.py")
import readwrite

def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Compute coexpression.")

    sys_args_parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="path to the 10X XeniumRanger output folder or proseg raw_results folder",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "--outcoexp",
        required=True,
        type=str,
        help="path to the output file for coexpression",
    )
    sys_args_parser.add_argument(
        "--outposrate",
        required=True,
        type=str,
        help="path to the output file for positive rate",
    )
    sys_args_parser.add_argument(
        "-m",
        required=True,
        type=str,
        help="method to compute coexpression",
    )
    sys_args_parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="target count for coexpression",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isdir(ret.i):
        raise RuntimeError(
            f"Error! 10X XeniumRanger output folder does not exist: {ret.i}"
        )
    os.makedirs(os.path.dirname(ret.outcoexp), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outposrate), exist_ok=True)

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # read counts
    adata = readwrite.read_xenium_sample(
        args.i,
        cells_as_circles=False,
        cells_boundaries=False,
        cells_boundaries_layers=False,
        nucleus_boundaries=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
        aligned_images=False,
        anndata=True,
    )

    # compute coexpression
    CC, X_downsampled, pos, pos_rate, mask = ce.coexpression(
        adata,
        target_count=args.c,
        method=args.m,
    )

    # save as parquet
    CC.to_parquet(args.outcoexp)
    pos_rate.to_frame().to_parquet(args.outposrate)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
