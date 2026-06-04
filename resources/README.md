# Resources

This directory contains supplementary material for Xenium data analysis outside
of the main pipeline. It includes tutorial notebooks and utility scripts for
image-related preprocessing workflows.

## Contents

- `notebooks/`: Jupyter notebooks for standalone analyses and tutorials.
- `scripts/`: Helper scripts for analysis preparation tasks.

## Notebooks

See `notebooks/README.md` for details.

- `notebooks/Xenium_Cell_Segmentation_Integration_with_QuPath.ipynb`: Tutorial
  notebook for enriching Xenium data with information from co-registered IF or
  H&E staining on the same slide.

## Scripts

See `scripts/README.md` for details.

- `scripts/README.md`: Usage notes for utility scripts.
- `scripts/cvt2ot/`: QuPath-based workflow for converting microscopy images,
  such as `.czi` or `.ndpi` H&E images, to pyramidal OME-TIFF files.
