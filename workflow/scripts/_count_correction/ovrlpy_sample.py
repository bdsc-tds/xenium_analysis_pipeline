import ovrlpy
import scipy
import sys
import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Compute ovrlpy correction.")
parser.add_argument(
    "--sample_transcripts_path",
    type=str,
    help="Path to the sample_signal_integrity file.",
)
parser.add_argument(
    "--out_file_signal_integrity",
    type=str,
    help="out_file_signal_integrity file to output.",
)
parser.add_argument(
    "--out_file_signal_strength",
    type=str,
    help="out_file_signal_strength file to output.",
)
parser.add_argument(
    "--out_file_transcript_info",
    type=str,
    help="out_file_transcript_info file to output.",
)
parser.add_argument(
    "--proseg_format",
    type="store_true",
    help="is the transcripts file in proseg raw output format.",
)
# parser.add_argument(
#     "--out_file_doublet_df",
#     type=str,
#     help="out_file_doublet_df file to output.",
# )

args = parser.parse_args()

# Access the arguments
sample_transcripts_path = args.sample_transcripts_path
out_file_signal_integrity = args.out_file_signal_integrity
out_file_signal_strength = args.out_file_signal_strength
out_file_transcript_info = args.out_file_transcript_info
proseg_format = args.proseg_format
# out_file_doublet_df = args.out_file_doublet_df


# load transcripts
if proseg_format:
    coordinate_df = (
        pd.read_csv(sample_transcripts_path, engine="pyarrow")
        .rename(
            columns={
                "assignment": "cell_id",
            }
        )
        .query("qv >= 20")
    )

    # remove dummy molecules
    coordinate_df = coordinate_df[
        ~coordinate_df["gene"].str.contains("|".join(["BLANK_", "UnassignedCodeword", "NegControl"]))
    ]
    # recode unassigned transcripts cell_id to UNASSIGNED
    coordinate_df["cell_id"] = (
        coordinate_df["cell_id"].astype(str).replace({str(coordinate_df["cell_id"].max()): "UNASSIGNED"})
    )

else:
    coordinate_df = (
        pd.read_parquet(sample_transcripts_path)
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

pd.DataFrame(signal_integrity).to_parquet(out_file_signal_integrity)
pd.DataFrame(signal_strength).to_parquet(out_file_signal_strength)
coordinate_df[["x_pixel", "y_pixel", "n_pixel", "z_delim"]].to_parquet(
    out_file_transcript_info
)

# doublet_df.to_parquet(out_file_doublet_df)
