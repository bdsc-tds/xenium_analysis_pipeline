import sys
import os
import argparse

import ovrlpy
import pandas as pd


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
        "--out_file_signal_integrity",
        type=str,
        required=True,
        help="out_file_signal_integrity file to output.",
    )
    parser.add_argument(
        "--out_file_signal_strength",
        type=str,
        required=True,
        help="out_file_signal_strength file to output.",
    )
    parser.add_argument(
        "--out_file_transcript_info",
        type=str,
        required=True,
        help="out_file_transcript_info file to output.",
    )
    parser.add_argument(
        "--proseg_format",
        action="store_true",
        help="is the transcripts file in proseg raw output format.",
    )
    # parser.add_argument(
    #     "--out_file_doublet_df",
    #     type=str,
    #     help="out_file_doublet_df file to output.",
    # )

    ret = parser.parse_args()
    if not os.path.isdir(ret.sample_transcripts_path):
        raise RuntimeError(
            f"Error! Input file does not exist: {ret.sample_transcripts_path}"
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
                "|".join(["BLANK_", "Codeword", "NegControl"])
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

    coordinate_df = coordinate_df.query("qv >= 20")  # remove low qv molecules
    coordinate_df["gene"] = coordinate_df["gene"].astype("category")

    # run ovrlpy
    signal_integrity, signal_strength, visualizer = ovrlpy.run(
        df=coordinate_df, cell_diameter=10, n_expected_celltypes=30
    )
    # doublet_df = ovrlpy.detect_doublets(
    #     signal_integrity, signal_strength, minimum_signal_strength=3, integrity_sigma=2
    # )

    # store results
    pd.DataFrame(signal_integrity).to_parquet(args.out_file_signal_integrity)
    pd.DataFrame(signal_strength).to_parquet(args.out_file_signal_strength)
    coordinate_df[["x_pixel", "y_pixel", "n_pixel", "z_delim"]].to_parquet(
        args.out_file_transcript_info
    )
    # doublet_df.to_parquet(args.out_file_doublet_df)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
