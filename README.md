# xenium_analysis_pipeline

## Reproducibility
This workflow is developed with reproducibility bearing in mind. Please refer to the following section for more details.

## Installation

### Snakemake
We use [Snakemake](https://snakemake.github.io) as the backend to this workflow. Thus, a conda environment for Snakemake must be created in the first step. We recommend [mamba](https://mamba.readthedocs.io/en/latest/index.html) as a replacement to conda for environment management.

Using `reproducibility/environment.yml`, we can create an environment for Snakemake:
```bash
# Alternative: `conda`
mamba env create -f reproducibility/environment.yml
```

### Singularity containers
We use multiple singularity containers for different methods and / or environments. To ensure reproducibility, please build these containers before executing the workflow.

#### R
The R version we use for this workflow is 4.4.1, and [renv](https://rstudio.github.io/renv/index.html) is used to track specific versions of packages. Please find files related to renv in `reproducibility/r/metadata`, and use `r.def` in `reproducibility/r` to build a corresponding container:
```bash
cd reproducibility/r
singularity build --fakeroot --force /path/to/the/built/container r.def
```

