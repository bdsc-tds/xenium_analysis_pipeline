#!/usr/bin/env bash
# Sourced by rule runSeggerNew (and the test sbatch) before calling `segger
# segment`. Fixes CUDA library discovery for the pip-installed RAPIDS/cupy stack
# when segger runs from a bare conda env (no container, no conda cudatoolkit):
#
#   1) cupy needs libcudart.so.12 etc. -> put the pip `nvidia-*-cu12` lib dirs on
#      LD_LIBRARY_PATH (they're not there by default).
#   2) cudf/numba dlopen the *unversioned* name (e.g. `libcudart.so`), which the
#      pip wheels don't ship -> create unversioned symlinks libX.so -> libX.so.NN.
#
# The env prefix is derived from segger's own location so this works for any env,
# including a Snakemake-rebuilt hashed-prefix conda env. Idempotent.
#
# NOTE: only the actual `segger segment` step needs this (it imports cupy/cudf);
# the segger_new2baysor converter does not.
_ep="$(dirname "$(dirname "$(command -v segger)")")"
_nvidia="$_ep"/lib/python*/site-packages/nvidia

export LD_LIBRARY_PATH="$(ls -d $_nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"

for _so in $_nvidia/*/lib/lib*.so.*; do
    [ -e "$_so" ] || continue
    ln -sf "$_so" "$(printf '%s' "$_so" | sed 's/\.so\..*/.so/')" 2>/dev/null || true
done
