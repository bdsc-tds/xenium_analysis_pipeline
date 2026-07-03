# Notebooks

This directory contains standalone notebooks for Xenium data analysis outside
of the main pipeline.

## Xenium Cell Segmentation Integration with QuPath

The notebook
`Xenium_Cell_Segmentation_Integration_with_QuPath.ipynb` describes how to enrich
Xenium data with information from co-registered IF or H&E staining on the same
slide.

## Requirements

The notebook requires an independent conda environment. Save the following
specification to a local file, for example `environment.yml`:

```yaml
name: img_registration
channels:
  - conda-forge
dependencies:
  - python=3.12
  - numpy=2.2.3
  - pandas=2.2.3
  - matplotlib=3.10.1
  - shapely=2.0.7
  - geopandas=1.0.1
  - jupyterlab=4.4.3
  - git
  - pip
  - pip:
      - git+https://github.com/scverse/spatialdata-io@4cf0f636a1ac78197e6f4fd92d4abafa5338ca61
      - spatialdata-plot==0.2.9
      - spatialdata==0.3.0
```

Create the environment with:

```bash
conda create -c conda-forge -f environment.yml
```

Activate the environment before starting JupyterLab:

```bash
conda activate img_registration
jupyter lab
```
