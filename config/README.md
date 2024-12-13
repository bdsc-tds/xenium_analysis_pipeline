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

- `_cell_type_annotation`: Disease specific configuration for cell type annotation using different approaches listed below. It should be specified as a direct child of the first layer (diseases).

  - `reference_based`: Single-cell RNA sequencing datasets are used as references for the annotation. Such reference datasets come either from a matched experiment (`matched_reference`), where samples are taken from the same tissue, or from an externel experiment (`external_reference`), where public datasets are used. Both types of references share the following specifications:

    - `path`: Absolute path to the reference datasets. An empty path (either `None` or `""`) indicates that the corresponding reference is unavailable, and thus the corresponding reference type won't be used in the pipeline for a specific disease. However, please note that at least one `path` of `matched_reference` and `external_reference` should be specified for annotating cell types.

    - `levels`: #TODO

- `_gene_panel_file`: The path or the name to the gene panel file. It should be specified as a direct child of the second layer (gene panels). Leaving blank implies that the gene panel is user defined and thus invariant to the different versions of 10X xeniumranger.

- `_qc`: Gene panel specific QC thresholds. It should be specified as a direct child of the second layer (gene panels). Only specify QC thresholds if the global ones are not suitable for specific gene panels. Please check the following section for global QC thresholds.

# config.yml

This file specifies the path to the experiment file and options for configuring tasks.

## Output path (key: `output_path`)

This section specifies the path to the output directory, which can either be relative with respect to the current working directory or be absolute.

## Singularity containers (key: `containers`)

This section specifies the paths to the prebuilt singularity containers, following the instructions in `README.md` in the root directory. The corresponding definition files can be found in `reproducibility` in the root directory.

## Experiments (key: `experiments`)

This section specifies the path, either absolute or relative, to the experiment file. If a relative path is provided, it must be relative to the `config.yml` file. Please refer to the above section (experiment.yml) for details.

## Reprocess (key: `reprocess`)

This section controls the step of reprocessing of each experiment. This step is useful when experiments, using default gene panels, were performed at different time points and thus processed by different versions of 10X xeniumranger. To eliminate the possible discrepancy due to renamed transcripts, this step reprocess every experiment with a specific version of xeniumranger. Whether to reprocess an experiment or not is determined by comparing the versions of 10X xeniumranger. The reprocessed results are then used as input to segementation.

Note that if customised gene panels are used, the results from 10X xeniumranger should be invariant to versions.

- `level`: The amount of discrepancy allowed between versions of 10X xeniumranger. Level 0 requires both versions to share the same major version, and level 1 requires both versions to share the same major and minor versions, and so on. Valid values are among [0, 1, 2, 3].

- `_threads`: The maximum number of threads to use.

- `_memory`: The maximum amount of memory (in GB) to use.

## Segmentation (key: `segmentation`)

This section specifies with key `methods` a number of methods that could be used for cell segmentation, including `10x`, `baysor`, `proseg`, and `segger`. Command line arguments for each of the methods are configured in a separate dictionary with the method name being the key. Removing any methods inside `methods` will rule them out from the workflow.

1. `10x`
   10X Xenium Ranger is used for segmentation. Users can specify some options of the method using the following keys (for the meaning for each option please refer to the official documentation):

   - `expansion-distance`: Either an integer or a list of integers. Each value is treated as an independent segmentation method, which is named as "10x\_{value}um".

   - `localcores`: The maximum number of threads to use. Either an integer or a list of integers. In the former case, it will be broadcasted into a list if the value of `expansion-distance` is a list.

   - `localmem`: The maximum amount of memory (in GB) to use. Similar to `localcores`.

   - `_other_options`: For options other than those above, in a form similar to: "--option_1 value_1 --option_2 value_2".

2. `baysor`
   [Baysor (v0.7.0+)](https://github.com/kharchenkolab/Baysor) is used for segmentation. Users can specify some options of the method using the following keys:

   - `_config`: The path, either absolute or relative, to the configuration file in TOML format for Baysor. If a relative path is given, it must be relative to the current working directory. By default the recommended configuration (stored in `workflow/configs/baysor_xenium.toml`) by developers of Baysor is used, assuming the current working directory is the root of this repo.

   - `_threads`: The maximum number of threads to use.

   - `_memory`: The maximum amount of memory (in GB) to use.

   - `_other_options`: For options not defined in the configuration, in a form similar to that mentioned above.

3. `proseg`
   [Proseg (v1.1.8+)](https://github.com/dcjones/proseg) is used for segmentation. Users can specify some options of the method using the following keys:

   - `_threads`: The maximum number of threads to use.

   - `_memory`: The maximum amount of memory (in GB) to use.

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
   - `max_counts`: `Inf` by default
   - `max_features`: `Inf` by default
   - `min_cells`: `5` by default

2. `dim_reduction`
   Parameters used in dimension reduction.

   - `n_dims`: Number of dimensions to keep. `50` by default.

3. `clustering`
   Parameters used in clustering.

   - `resolution`: Parameter for `Seurat::FindClusters`. `0.8` by default.

In addition, there are some reserved keys, specified as follows:

1. `_future_globals_maxSize`: A number (`1` by default) with unit GB is specified here to deal with the [issue](https://github.com/satijalab/seurat/issues/1845) of "Global size exceeds maximum allowed size" when running Seurat (see the solution [here](https://satijalab.org/seurat/archive/v3.0/future_vignette.html)).

## Cell type annotation (key: `cell_type_annotation`)

This section specifies with key `approaches` some approaches to use for cell type annotation. For now only reference based approach (specified with key `reference_based`) is available.

More detailed configuration for each approach is configured in a separate dictionary with the approach name being the key. Removing any approaches inside `approaches` will rule them out from the workflow.

1. `reference_based`
   Use reference-based methods for cell type annotation, where disease-specific single-cell RNA sequencing datasets are used as references (please refer to the details in section `experiments.yml` above).

   With key `methods`, this section specifies a number of methods that could be used for cell type annotation, including `rctd`, `singler`, `seurat`, `tangram`, and `xgboost`. Command line arguments for each of the methods are configured in a separate dictionary with the method name being the key. Removing any methods inside `methods` will rule them out from the workflow.

   1.1 `rctd`

   1.2 `singler`

   1.3 `seurat`

   1.4 `tangram`

   1.5 `xgboost`
