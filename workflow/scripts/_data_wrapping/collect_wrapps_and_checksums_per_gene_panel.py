import sys

old_stdout = sys.stdout
old_stderr = sys.stderr
_log = open(snakemake.log[0], "w", encoding="utf-8")
sys.stdout = _log
sys.stderr = _log

from pathlib import Path
import subprocess
import pandas as pd

# Combine wrapped samples into a single tarball per gene panel
print(
    "Combining wrapped samples into a single tarball per gene panel...\nSamples to be combined:"
)
print("\n".join(snakemake.input.wrap))

cmd = [
    "tar",
    "-czf",
    snakemake.output.wrap,
    "-C",
    snakemake.params.work_dir,
    *[Path(i).name for i in snakemake.input.wrap],
]

print(f"Running command: {' '.join(cmd)}")

subprocess.run(
    cmd,
    check=True,
)

# Combine checksum files into a single checksum file per gene panel
print(
    "Combining checksum files into a single checksum file per gene panel...\nChecksum files to be combined:"
)
print("\n".join(snakemake.input.checksum))

comb_cs = None
for i in snakemake.input.checksum:
    df = pd.read_csv(
        i,
        sep=",",
        header=0,
        index_col=False,
    )
    if comb_cs is None:
        comb_cs = df
    else:
        comb_cs = pd.concat(
            [comb_cs, df],
            ignore_index=True,
        )

print(comb_cs)

comb_cs.to_csv(
    snakemake.output.checksum,
    sep=",",
    header=True,
    index=False,
)

_log.close()
sys.stdout = old_stdout
sys.stderr = old_stderr
