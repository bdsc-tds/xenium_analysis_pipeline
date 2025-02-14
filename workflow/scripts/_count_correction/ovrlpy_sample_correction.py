import argparse
import sys
import os

import pandas as pd
import anndata as ad
from ..utils import readwrite


def transcripts_to_count_matrix(
    transcripts, cell_column="cell_id", feature_column="feature_name", qv_treshold=20
):
    transcripts = transcripts.query(
        f"is_gene &(qv >= {qv_treshold}) & ({cell_column} != 'UNASSIGNED')"
    )
    cm = transcripts.pivot_table(
        index=cell_column, columns=feature_column, aggfunc="size", fill_value=0
    )
    return cm


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Compute ovrlpy correction.")
    parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument(
        "--sample_transcripts_path",
        type=str,
        required=True,
        help="Path to the sample_signal_integrity file.",
    )
    parser.add_argument(
        "--sample_signal_integrity",
        type=str,
        required=True,
        help="Path to the sample_signal_integrity file.",
    )
    parser.add_argument(
        "--sample_transcript_info",
        type=str,
        required=True,
        help="sample_transcript_info file to output.",
    )
    parser.add_argument(
        "--out_file_corrected_counts",
        type=str,
        required=True,
        help="out_file_corrected_counts file to output.",
    )
    parser.add_argument(
        "--out_file_cells_mean_integrity",
        type=str,
        required=True,
        help="out_file_corrected_counts file to output.",
    )
    parser.add_argument(
        "--signal_integrity_threshold",
        type=float,
        required=True,
        help="signal_integrity_threshold parameter (threshold below which a pixel is low quality).",
    )
    parser.add_argument(
        "--proseg_format",
        action="store_true",
        help="is the transcripts file in proseg raw output format.",
    )

    ret = parser.parse_args()
    if not os.path.isdir(ret.sample_transcripts_path):
        raise RuntimeError(
            f"Error! Input file does not exist: {ret.sample_transcripts_path}"
        )
    if not os.path.isdir(ret.sample_signal_integrity):
        raise RuntimeError(
            f"Error! Input file does not exist: {ret.sample_signal_integrity}"
        )
    if not os.path.isdir(ret.sample_transcript_info):
        raise RuntimeError(
            f"Error! Input file does not exist: {ret.sample_transcript_info}"
        )
    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # load transcripts
    if args.proseg_format:
        coordinate_df = (
            pd.read_csv(args.sample_transcripts_path, engine="pyarrow")
            .rename(
                columns={
                    "assignment": "cell_id",
                }
            )
            .query("qv >= 20")
        )

        # remove dummy molecules
        coordinate_df = coordinate_df[
            ~coordinate_df["gene"].str.contains(
                "|".join(["BLANK_", "UnassignedCodeword", "NegControl"])
            )
        ]
        # recode unassigned transcripts cell_id to UNASSIGNED
        coordinate_df["cell_id"] = (
            coordinate_df["cell_id"]
            .astype(str)
            .replace({str(coordinate_df["cell_id"].max()): "UNASSIGNED"})
        )
    else:
        coordinate_df = (
            pd.read_parquet(args.sample_transcripts_path)
            .rename(
                columns={
                    "x_location": "x",
                    "y_location": "y",
                    "z_location": "z",
                    "feature_name": "gene",
                }
            )
            .query("is_gene")  # remove dummy molecules
        )
    transcript_info = pd.read_parquet(args.sample_transcript_info)
    coordinate_df = coordinate_df.join(transcript_info)

    # load ovrlpy
    signal_integrity = pd.read_parquet(args.sample_signal_integrity).values
    # x and y are transposed in ovrlpy output. sanity check x and y must be transposed
    assert (coordinate_df.y_pixel.max() + 1) == signal_integrity.shape[0]
    assert (coordinate_df.x_pixel.max() + 1) == signal_integrity.shape[1]
    coordinate_df["signal_integrity"] = signal_integrity[
        coordinate_df.y_pixel, coordinate_df.x_pixel
    ]
    # filter transcripts based on ovrlpy signal integrity
    coordinate_df_filtered = coordinate_df[
        coordinate_df.signal_integrity > args.signal_integrity_threshold
    ]

    # compute mean cell integrity QC metric from unfiltered transcripts
    cell_mean_integrity = (
        coordinate_df.query("cell_id != 'UNASSIGNED'")
        .groupby("cell_id")["signal_integrity"]
        .mean()
    )
    # cell_is_singlet = cell_mean_integrity > args.signal_integrity_threshold

    # create filtered count matrix
    corrected_counts = transcripts_to_count_matrix(
        coordinate_df_filtered, feature_column="gene"
    )

    # store results
    adata_out = ad.AnnData(corrected_counts)
    readwrite.write_10X_h5(adata_out, args.out_file_corrected_counts)
    if args.out_file_cells_mean_integrity != "":
        cell_mean_integrity.to_frame().to_parquet(args.out_file_cells_mean_integrity)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
