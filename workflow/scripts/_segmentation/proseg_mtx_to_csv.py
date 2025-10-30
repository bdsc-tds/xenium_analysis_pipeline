import sys
import pandas as pd
from pathlib import Path
from scipy.io import mmread


def main():
    p = Path(sys.argv[1])

    mtx_file = p / "expected-counts.mtx.gz"
    gene_file = p / "gene-metadata.csv.gz"
    out_file = p / "expected-counts.csv.gz"

    print("Reading matrix")
    X = mmread(mtx_file).toarray()

    print("Reading gene info")
    genes = pd.read_csv(gene_file)["gene"].values
    df = pd.DataFrame(X, columns=genes)

    print("writing csv")
    df.to_csv(out_file, compression="gzip", index=False)


if __name__ == "__main__":
    main()
