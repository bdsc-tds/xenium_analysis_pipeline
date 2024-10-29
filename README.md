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

We have developed a bash script, `run.sh`, to make the execution easy for users. Before execution, users need to configure a few entries in it under _USER SETUP_ section.

- `MODULES`: Modules to be loaded prior to execution on clusters.

- `CONDA_BIN`: The name of or path to either `mamba` or `conda`.

- `ENV_NAME`: The name of or path to the conda environment (by default `xenium_analysis_pipeline`).

- `LOCAL_PROFILE` and `CLUSTER_PROFILE`:
  The execution of this workflow is controlled by profiles. Please refer to [the Snakemake manuel](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for the details.

  We have provided two examples of profiles under `profiles`. One is for local execution, which locates in `profiles/local`; the other is for cluster execution, specifically slurm, which locates in `profiles/slurm`. Users can edit these profiles according to their specific needs.

  Besides, users can define their own profiles for execution, e.g., when they use a cluster other than slurm. They only need to specify in proper places the paths to their customised profiles.

- `SINGULARITY_BIND_DIRS`: An array of directories to bind to containers. Each element should be in the following form: _LOCAL_DIR:SINGULARITY_DIR_. Inexistent local directories will be filtered out.

## Execution

The workflow should be executed from the root directory of this repo. To get a self-explanatory help message, type

```bash
# the current working directory is the root of this repo
./run.sh --help
```

which prints

```
Usage: [ -m | --mode MODE ] [ -c | --core CORE ] [ -n | --dry-run ] [ --dag OUTPUT ] [ --unlock ] [ -v | --verbose ] [ -h | --help ]
        -m,--mode MODE: the pipeline will be run on 'local' (default) or on 'cluster'.
        -c,--core CORE: the number of cores to be used when -m,--mode is unset or 'local' (default: 1); ignored when -m,--mode is 'cluster'.
        -n,--dry-run: dry run.
        --dag OUTPUT: draw dag and save to OUTPUT.pdf.
        --unlock: unlock the working directory.
        -v,--verbose: print more information.
        -h,--help: print this message.
```

## Solutions to known problems

1. Snakemake fails to create conda environments due to the lack of writing permission of `/tmp/conda`.

   Such an issue occurs most often on devices with multiple users, such as HPCs and servers. In this case, the reason is simply that `/tmp/conda` is possessed by other users. A possible solution is to firstly create a folder in `/tmp`, such as `/tmp/your_id`, and then add `/tmp/your_id:/tmp` to `SINGULARITY_BIND_DIRS` inside `run.sh`. After this you can rerun the command, and those environments should be correctly created.

   Additionally, for HPCs, users might need to remove `/tmp/your_id:/tmp` from `SINGULARITY_BIND_DIRS` after creating environments as directory `/tmp/your_id` is not likely present in compute nodes.

2. For those steps involving 10X xeniumranger, sometimes I get the folowing error: "PermissionError: [Errno 13] Permission denied".

   10X xeniumranger copies files from raw data during processing. This error could be because the user, as the owner of the raw data, deprives him-/herself of write permission to it. When 10X xeniumranger conducts copy operation, it also copies the modes of files, and hence this error when it needs to write to the copied files. Although it is a safe behaviour to prevent from accidental change of the raw data, users have to have write permission to the raw data when they are also the owner.
