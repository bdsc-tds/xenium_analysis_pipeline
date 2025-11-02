import argparse
import os
import sys
import pandas as pd
from scipy.io import mmread


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert ProSeq matrix file to CSV format."
    )
    parser.add_argument(
        "-l",
        type=str,
        help="Path to log file.",
    )
    parser.add_argument(
        "--counts",
        type=str,
        required=True,
        help="Path to the input counts file.",
    )
    parser.add_argument(
        "--expected_counts",
        type=str,
        required=True,
        help="Path to the input expected counts file.",
    )
    parser.add_argument(
        "--gene",
        type=str,
        required=True,
        help="path to the input gene metadata file",
    )
    parser.add_argument(
        "--out_counts",
        type=str,
        required=True,
        help="Path to the output counts CSV file.",
    )
    parser.add_argument(
        "--out_expected_counts",
        required=True,
        help="Path to the output expected counts CSV file.",
    )

    ret = parser.parse_args()
    if not os.path.isfile(ret.counts):
        raise RuntimeError(
            f"Error! The provided counts file does not exist: {ret.counts}"
        )
    if not os.path.isfile(ret.expected_counts):
        raise RuntimeError(
            f"Error! The provided expected counts file does not exist: {ret.expected_counts}"
        )

    os.makedirs(os.path.dirname(ret.out_counts), exist_ok=True)
    os.makedirs(os.path.dirname(ret.out_expected_counts), exist_ok=True)

    return ret


def _convert(
    gene_info,
    mtx_file: str,
    out_file: str,
) -> None:
    print("Reading matrix")
    X = mmread(mtx_file).toarray()

    df = pd.DataFrame(X, columns=gene_info)

    print("writing csv")
    df.to_csv(out_file, compression="gzip", index=False)


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    print("Reading gene info")
    genes = pd.read_csv(args.gene)["gene"].values

    _convert(
        genes,
        args.counts,
        args.out_counts,
    )
    _convert(
        genes,
        args.expected_counts,
        args.out_expected_counts,
    )

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
