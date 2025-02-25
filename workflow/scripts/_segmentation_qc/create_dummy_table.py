"""
This script is used to generate a dummy table.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Generate a dummy table.")

    sys_args_parser.add_argument(
        "--segmentation_id",
        required=True,
        type=str,
        help="segmentation id",
    )
    sys_args_parser.add_argument(
        "--segmentation_method",
        required=True,
        type=str,
        help="segmentation method",
    )
    sys_args_parser.add_argument(
        "--sample_id",
        required=True,
        type=str,
        help="sample id",
    )
    sys_args_parser.add_argument(
        "--condition",
        required=True,
        type=str,
        help="condition",
    )
    sys_args_parser.add_argument(
        "--gene_panel",
        required=True,
        type=str,
        help="gene panel",
    )
    sys_args_parser.add_argument(
        "--donor",
        required=True,
        type=str,
        help="donor",
    )
    sys_args_parser.add_argument(
        "--sample",
        required=True,
        type=str,
        help="sample",
    )
    sys_args_parser.add_argument(
        "--normalisation_id",
        required=True,
        type=str,
        help="normalisation id",
    )
    sys_args_parser.add_argument(
        "--annotation_id",
        type=str,
        help="annotation id",
    )
    sys_args_parser.add_argument(
        "--reference_name",
        required=True,
        type=str,
        help="reference name",
    )
    sys_args_parser.add_argument(
        "--reference_level",
        type=str,
        help="reference level",
    )
    sys_args_parser.add_argument(
        "--na_values",
        action="store_true",
        help="whether to use nan for missing values (0 as alternative)",
    )
    sys_args_parser.add_argument(
        "--out",
        required=True,
        type=str,
        help="path to the output file",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )

    ret = sys_args_parser.parse_args()
    os.makedirs(os.path.dirname(ret.out), exist_ok=True)

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    missing_value = np.nan if args.na_values else 0

    pd.DataFrame(
        {
            "segmentation_id": args.segmentation_id,
            "segmentation_method": args.segmentation_method,
            "sample_id": args.sample_id,
            "condition": args.condition,
            "gene_panel": args.gene_panel,
            "donor": args.donor,
            "sample": args.sample,
            "normalisation_id": args.normalisation_id,
            "annotation_id": args.annotation_id,
            "reference_name": args.reference_name,
            "reference_level": args.reference_level,
            "mean_ncount": missing_value,
            "median_ncount": missing_value,
            "mean_nfeature": missing_value,
            "median_nfeature": missing_value,
            "ncell": missing_value,
            "reject": missing_value,
            "singlet": missing_value,
            "doublet_certain": missing_value,
            "doublet_uncertain": missing_value,
            "unannotated": missing_value,
            "reject_prop": missing_value,
            "singlet_prop": missing_value,
            "doublet_certain_prop": missing_value,
            "doublet_uncertain_prop": missing_value,
            "unannotated_prop": missing_value,
        },
        index=[0],
    ).to_parquet(
        args.out,
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
