#!/bin/bash

#SBATCH --job-name=qupath_conv
#SBATCH --output=logs/conv_%a.out
#SBATCH --error=logs/conv_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --time=04:00:00

# Change this to reflect the actual number of images to be converted.
#SBATCH --array=1-1

# User setup
# Directory to images as input
INPUT_DIR=
# Directory to the converted images as output
OUTPUT_DIR=
# Pattern of files to be converted, e.g., "*.czi", "*.ndpi", etc.
FILE_PATTERN="*"
# Path to QuPath binary
QUPATH_BIN=QuPath

# DO NOT EDIT BELOW
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: INPUT_DIR does not exist or is not a directory: $INPUT_DIR" >&2
  exit 1
fi

if [[ "$QUPATH_BIN" == */* ]]; then
  if [[ ! -e "$QUPATH_BIN" ]]; then
    echo "Error: QUPATH_BIN does not exist: $QUPATH_BIN" >&2
    exit 1
  fi
else
  if ! command -v "$QUPATH_BIN" >/dev/null 2>&1; then
    echo "Error: QUPATH_BIN was not found in PATH: $QUPATH_BIN" >&2
    exit 1
  fi
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get the list of files and pick the one for THIS task
# Files and store them in an array
IFS=$'\n' FILE_LIST=( $(find "${INPUT_DIR}" -maxdepth 1 -name "${FILE_PATTERN}" | sort) )
# Note: Arrays are 0-indexed, SLURM tasks usually start at 1
CURRENT_FILE="${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

echo "Task $SLURM_ARRAY_TASK_ID processing ${CURRENT_FILE}"

# Run QuPath
"${QUPATH_BIN}" -Djava.awt.headless=true script \
  -a output \
  -a "${OUTPUT_DIR}" \
  -i "${CURRENT_FILE}" \
  cvt2ot.groovy

echo "Done!"
