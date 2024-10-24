# experiments.yml

This file lists experiments to be analysed.

Experiment data are expected to be organised in the following manner on disk, where each layer is a folder:

- Layer 1: diseases
- Layer 2: gene panels
- Layer 3: donors
- Layer 4: samples

There is no specific naming convention for each layer, except that they should not start with underscores, as character strings with such patterns are reserved for special usage by this workflow.

With the data organised, they should then be specified in `experiments.yml`. The first three layers, namely diseases, gene panels and donors, should be listed as nested keys, and the final layer (samples) is a list of values corresponding to the third layer (donors). Users must ensure that these keys and values are character strings.

Some reserved keys are listed as follows:

- `_base_path`: The parent directory to the first layer (diseases) of the organised data. It should be specified as siblings of the first layer (diseases).

- `_gene_panel_file`: The path or the name to the gene panel file. It should be specified as the direct child of the second layer (gene panels). Leaving blank implies that the gene panel is user defined and thus invariant to the different versions of 10X xeniumranger.

- `_qc`: Gene panel specific QC thresholds. Only specify QC thresholds if the global ones are not suitable for specific gene panels. Check the folling section for global QC thresholds.

# config.yml

This file specifies the path to the experiment file and options for configuring tasks.

## Output path (key: `output_path`)

This section specifies the path to the output directory, which can either be relative with respect to the current working directory or be absolute.

## Singularity containers (key: `containers`)

This section specifies the paths to the prebuilt singularity containers, following the instructions in `README.md` in the root directory. The corresponding definition files can be found in `reproducibility` in the root directory.

## Experiments (key: `experiments`)

This section specifies the path, either absolute or relative, to the experiment file. If a relative path is provided, it must be relative to the `config.yml` file. Please refer to the above section (experiment.yml) for details.

## Segmentation (key: `segmentation`)

This section specifies with key `methods` a number of methods that could be used for cell segmentation, including `10x`, `baysor`, `proseg`, and `segger`. Command line arguments for each of the methods are configured in a separate dictionary with the method name being the key. Removing any methods inside `methods` will rule them out from the workflow.

1. `10x`

10X Xenium Ranger is used for segmentation. Users can specify some options of the method using the following keys (for the meaning for each option please refer to the official documentation):

- `expansion-distance`: Either an integer or a list of integers. Each value is treated as an independent segmentation method, which is named as "10x\_{value}um".

- `localcores`: The maximum number of threads to use. Either an integer or a list of integers. In the former case, it will be broadcasted into a list if the value of `expansion-distance` is a list.

- `localmem`: The maximum amount of memory (in GB) to use. Similar to `localcores`.

- `_other_options`: For options other than those above, in a form similar to: "--option_1 value_1 --option_2 value_2".

2. `baysor`

[Baysor](https://github.com/kharchenkolab/Baysor) is used for segmentation. Users can specify some options of the method using the following keys:

- `_config`: The path, either absolute or relative, to the configuration file in TOML format for Baysor. If a relative path is given, it must be relative to the current working directory. By default the recommended configuration (stored in `workflow/configs/baysor_xenium.toml`) by developers of Baysor is used, assuming the current working directory is the root of this repo.

- `_other_options`: For options not defined in the configuration, in a form similar to that mentioned above.

In addition, after acquiring results from Baysor, Proseg, and Segger, an extra step is performed where those results are normalised into the same format for the sake of downstream analysis. The corresponding arguments are specified under `_normalisation`, including:

- `localcores`: The maximum number of threads to use.

- `localmem`: The maximum amount of memory (in GB) to use.

## Standard Seurat analysis (key: `standard_seurat_analysis`)

This section specifies parameters in the standard Seurat analysis workflow.

1. `qc`

Global QC thresholds. If any of the following QC thresholds are not suitable for certain gene panels, users should specify them accordingly in `experiments.yml`.

- `min_counts`: `10` by default.
- `min_features`: `5` by default
- `max_counts`: `.inf` by default
- `max_features`: `.inf` by default
- `min_cells`: `1` by default
