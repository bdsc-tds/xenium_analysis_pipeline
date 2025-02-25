import sys

import numpy as np
import pandas as pd

old_stdout = sys.stdout
old_stderr = sys.stderr
_log = open(snakemake.log[0], "w", encoding="utf-8")
sys.stdout = _log
sys.stderr = _log


def per_row_func(row: pd.Series) -> pd.Series:
    if np.isnan(row.ncell):
        return row

    if row.ncell > 0:
        return row.fillna(
            value=0,
        )

    return row


(
    pd.concat(
        [
            pd.read_parquet(
                i,
            )
            for i in snakemake.input
        ],
    )
    .reset_index(
        drop=True,
    )
    .apply(
        per_row_func,
        axis=1,
    )
    .to_parquet(
        snakemake.output[0],
    )
)


_log.close()
sys.stdout = old_stdout
sys.stderr = old_stderr
