# Utility Scripts

This directory contains helper scripts related to Xenium data analysis outside of the pipeline.

## Convert Images to OME-TIFF with QuPath

`cvt2ot` contains scripts converting image files, such as H&E images of `.czi` or `.ndpi` format,
to pyramidal OME-TIFF files using QuPath in headless mode.

Files:

- `cvt2ot/run.sh`: SLURM array job wrapper.
- `cvt2ot/cvt2ot.groovy`: QuPath script that exports each image or image
  series to pyramidal `.ome.tif`.

### Requirements

- QuPath command line executable available either in `PATH` or as an explicit
  path in `QUPATH_BIN`.
- A SLURM cluster if submitting `run.sh` with `sbatch`.
- Input image files readable by QuPath/Bio-Formats.

### Configure `run.sh`

Edit `resources/scripts/cvt2ot/run.sh` before submitting the job:

```bash
# Directory to images as input
INPUT_DIR=/path/to/input_images

# Directory to the converted images as output
OUTPUT_DIR=/path/to/output_ome_tiff

# Pattern of files to be converted, e.g., "*.czi", "*.ndpi", etc.
FILE_PATTERN="*"

# Path to QuPath binary
QUPATH_BIN=QuPath
```

If QuPath is not available in `PATH`, set `QUPATH_BIN` to the executable path,
for example:

```bash
QUPATH_BIN=/path/to/QuPath/bin/QuPath
```

The script checks that `INPUT_DIR` exists and that `QUPATH_BIN` is available
before creating output directories or starting conversion.

### Set the SLURM Array Size

`run.sh` processes one input file per SLURM array task. Update this line to
match the number of files selected by `FILE_PATTERN`:

```bash
#SBATCH --array=1-1
```

For example, if there are 12 `.czi` files:

```bash
#SBATCH --array=1-12
```

You can count matching files with:

```bash
find /path/to/input_images -maxdepth 1 -name "*.czi" | sort | wc -l
```

### Submit the Job

Submit the job from the `cvt2ot` directory because `run.sh` references
`cvt2ot.groovy` by relative path:

```bash
cd resources/scripts/cvt2ot
sbatch run.sh
```

SLURM logs are written to:

```text
resources/scripts/cvt2ot/logs/conv_<array_task_id>.out
resources/scripts/cvt2ot/logs/conv_<array_task_id>.err
```

### Output Naming

For each input image, `cvt2ot.groovy` writes pyramidal OME-TIFF files to
`OUTPUT_DIR`.

If an image has a single large series, the output name is:

```text
<input_basename>.ome.tif
```

If an image contains multiple series or named regions, the script appends a
sanitized series or metadata suffix:

```text
<input_basename>_<series_or_region_name>.ome.tif
```

Small utility images, such as labels or thumbnails, are skipped when either
dimension is below 1000 pixels.

### Export Settings

The Groovy script exports with:

- OME pyramidal TIFF format.
- Downsamples: `1, 2, 4, 8, 16, 32`.
- Tile size: `1024`.
- Compression: `ZLIB`.
- Parallel writing enabled.

### Troubleshooting

- `Error: INPUT_DIR does not exist or is not a directory`: set `INPUT_DIR` to an
  existing image directory.
- `Error: QUPATH_BIN was not found in PATH`: load the QuPath module or set
  `QUPATH_BIN` to the full QuPath executable path.
- No files are processed: check `FILE_PATTERN` and make sure the SLURM array
  size is not larger than the number of matching files.
- `cvt2ot.groovy` cannot be found: submit the job from
  `resources/scripts/cvt2ot`, or update `run.sh` to pass the full path to the
  Groovy script.
