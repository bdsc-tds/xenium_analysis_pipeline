import sys

import pandas as pd

old_stdout = sys.stdout
old_stderr = sys.stderr
_log = open(snakemake.log[0], "w", encoding="utf-8")
sys.stdout = _log
sys.stderr = _log


info: pd.DataFrame = (
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
    .fillna(
        value=0,
    )
)

info.to_parquet(
    snakemake.output[0],
)


_log.close()
sys.stdout = old_stdout
sys.stderr = old_stderr
