# xenium_analysis_pipeline

## Reproducibility

This workflow is developed with reproducibility bearing in mind. Please refer to the following section for more details.

## Installation

### Snakemake

We use [Snakemake (v8+)](https://snakemake.github.io) as the backend to this workflow. Thus, a conda environment for Snakemake must be created in the first step. We recommend [mamba](https://mamba.readthedocs.io/en/latest/index.html) as a replacement to conda for environment management.

Using `reproducibility/environment.yml`, we can create an environment for Snakemake:

```bash
# the current working directory is the root of this repo
# Alternative: `conda`
mamba env create -f reproducibility/environment.yml
```

### Singularity containers

We use multiple singularity containers for different methods and / or environments. To ensure reproducibility, please build these containers before executing the workflow.

#### R

The R version we use for this workflow is 4.4.1, and [renv](https://rstudio.github.io/renv/index.html) is used to track specific versions of packages. Please find files related to renv in `reproducibility/r/metadata`, and use `r.def` in `reproducibility/r` to build the corresponding container:

```bash
# the current working directory is the root of this repo
cd reproducibility/r
singularity build --fakeroot --force /path/to/the/built/container r.def
```

#### 10X Xenium Ranger

The 10X Xenium Ranger version we use here is 3.0.1 (Sep. 19th, 2024). A link is used to download the software from the [10X website](https://www.10xgenomics.com/support/software/xenium-ranger/downloads). Since 10X regularly updates this link, users should replace it with the most recent one if the container fails to be built:

```bash
# the current working directory is the root of this repo
singularity build --fakeroot --force /path/to/the/built/container reproducibility/10x.def
```

#### Baysor

The [Baysor](https://github.com/kharchenkolab/Baysor) version we use here is 0.7.0.

```bash
# the current working directory is the root of this repo
singularity build --fakeroot --force /path/to/the/built/container reproducibility/baysor.def
```

## Configuration

### Configuration for the workflow

Please edit `config/config.yml` for the configuration of the workflow. A detailed guideline can be found in `config/README.md`.

### Configuration for execution

The execution of this workflow is controlled by profiles. Please refer to [the Snakemake manuel](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for the details.

We have provided two examples of profiles under `profiles`. One is for local execution, which locates in `profiles/local`; the other is for cluster execution, specifically slurm, which locates in `profiles/slurm`. Users can edit these profiles according to their specific needs.

Besides, users can define their own profiles for execution. They only need to specify the paths to their customised profiles in `run.sh`. See the section below for more details.
