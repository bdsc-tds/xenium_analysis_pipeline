import sys

old_stdout = sys.stdout
old_stderr = sys.stderr

_log = open(snakemake.log[0], "w", encoding="utf-8")

sys.stdout = _log
sys.stderr = _log


import json

with open(snakemake.output[0], "w", encoding="utf-8") as f:
    json.dump(
        snakemake.params[0],
        f,
    )


_log.close()
sys.stdout = old_stdout
sys.stderr = old_stderr
