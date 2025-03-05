import sys
import os
import argparse

import pandas as pd

sys.path.append(
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
        ),
    ),
)
from utils import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Compute ovrlpy correction.")
    parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument(
        "-i",
        type=str,
        required=True,
        help="Path to the input data directory.",
    )
    parser.add_argument(
        "-o",
        type=str,
        required=True,
        help="Path to the output file.",
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        required=True,
        help="meta data - sample_id",
    )
    parser.add_argument(
        "--segmentation_id",
        type=str,
        required=True,
        help="meta data - segmentation_id",
    )
    parser.add_argument(
        "--condition",
        type=str,
        required=True,
        help="meta data - condition",
    )
    parser.add_argument(
        "--gene_panel",
        type=str,
        required=True,
        help="meta data - gene_panel",
    )
    parser.add_argument(
        "--donor",
        type=str,
        required=True,
        help="meta data - donor",
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="meta data - sample",
    )
    parser.add_argument(
        "--segmentation_method",
        type=str,
        required=True,
        help="meta data - segmentation_method",
    )
    parser.add_argument(
        "--in_mapping",
        type=str,
        help="Path to the input mapping file.",
    )
    parser.add_argument(
        "--cell_id_col_name",
        type=str,
        help="cell id column name",
    )
    parser.add_argument(
        "--use_mapping",
        action="store_true",
        help="whether to use mapping, which uses a subset of cells",
    )

    ret = parser.parse_args()
    if not os.path.isdir(ret.i):
        raise RuntimeError(f"Error! Input directory does not exist: {ret.i}")

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

    data: readwrite.ad.AnnData = readwrite.read_xenium_sample(
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

    data = readwrite.add_metadata2ad(
        data,
        sample_id=args.sample_id,
        segmentation_id=args.segmentation_id,
        condition=args.condition,
        gene_panel=args.gene_panel,
        donor=args.donor,
        sample=args.sample,
        segmentation_method=args.segmentation_method,
    )

    if (
        args.in_mapping is not None
        and args.cell_id_col_name is not None
        and args.use_mapping
    ):
        assert os.path.isfile(
            args.in_mapping,
        )

        mapping = pd.read_parquet(
            args.in_mapping,
        )

        data = data[
            mapping[args.cell_id_col_name],
            :,
        ]

    data.write(
        args.o,
        compression="gzip",
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
