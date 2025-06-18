This directory contains several notebooks related to Xenium data analysis outside of the pipeline.

# Xenium cell segmentation integration with QuPath

This analysis enriches the Xenium data with insights from co-registered IF or H&E staining on the same slide. It has been described in detail in the notebook named `Xenium_Cell_Segmentation_Integration_with_QuPath.ipynb`.

The execution of the code blocks in the tutorial requires an independent virtual environment. Please save the following specifications for the conda environment to a local file named, for instance, `environment.yml`.

```
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

Then create the environment with `conda create -c conda-forge -f environment.yml`.
