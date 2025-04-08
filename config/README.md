# experiments.yml

This file lists experiments to be analysed.

Experiment data are expected to be organised in the following manner on disk, where each layer is a folder:

- Layer 1: conditions
- Layer 2: gene panels
- Layer 3: donors
- Layer 4: samples

There is no specific naming convention for each layer, except that they should not contain pluses (+) or start with underscores, as character strings with such patterns are reserved for special usage by this workflow.

With the data organised, they should then be specified in `experiments.yml`. The first three layers, namely conditions, gene panels and donors, should be listed as nested keys, and the final layer (samples) is a list of values corresponding to the third layer (donors). Users must ensure that these keys and values are character strings.

Some reserved keys are listed as follows:

- `_base_path`: The parent directory to the first layer (conditions) of the organised data. It should be specified as siblings of the first layer (conditions).

- `_cell_type_annotation`: Condition specific configuration for cell type annotation using different approaches listed below. It should be specified as a direct child of the first layer (conditions).

  - `reference_based`: Annotated single-cell RNA sequencing datasets are used as references. Such reference datasets come either from a matched experiment, where samples are taken from the same tissue, or from an externel experiment, where public datasets are used. Users can specify multiple references with unique, user-defined names (not starting with underscores), which share the following specifications:

    - `path`: Absolute path to the reference dataset. An empty path (either `None` or `""`) indicates that the corresponding reference is unavailable, and thus the corresponding reference type won't be used in the pipeline for a specific condition. Note that it is required for the reference dataset to have a column named "nCount", representing the number of counts per cell.

    - `levels`: A list of level names in the reference dataset. A valid reference dataset also requires a valid list of levels.

    - `cell_min_instance`: The minimum number of cells required per cell type (by default `25`). Can be decreased for smaller references.

- `_gene_panel_file`: The path or the name to the gene panel file. It should be specified as a direct child of the second layer (gene panels). Leaving blank implies that the gene panel is user defined and thus invariant to the different versions of 10X xeniumranger.

- `_qc`: Gene panel specific QC thresholds. It should be specified as a direct child of the second layer (gene panels). Only specify QC thresholds if the global ones are not suitable for specific gene panels. Please check the following section for global QC thresholds.

- `_target_counts`: gene panel specific target counts for computing coexpression values. It should be specified as a direct child of the second layer (gene panels). Only specify target values if the global ones are not suitable for specific gene panels. Please check the following section for global target counts.

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

This section specifies with key `methods` a number of methods that could be used for cell segmentation, including `10x_mm`, `10x`, `baysor`, `proseg_expected`, `proseg_mode`, and `segger`. Specifying either or both of `proseg_expected` and `proseg_mode` will only run `proseg` once. The only difference is that which read counts are used in the downstream analyses: for `proseg_expected`, the expected read counts are used, as recommended by the author of `proseg`; for `proseg_mode`, the read counts are computed by aggregating transcripts assigned to cells based on the mode of their posterior distribution.

Command line arguments for each of the methods are configured in a separate dictionary with the method name being the key. Removing any methods inside `methods` will rule them out from the workflow.

1. `10x_mm`
   10X Xenium Ranger's multimodal cell segmentation algorithm is used, where options `--boundary-stain` and `--interior-stain` are selected automatically based on present stain data. This method should only be activated when the boundary and interior stain signals are available in the data.

   Users can specify some options of the method using the following keys (for the meaning for each option please refer to the official documentation):

   - `expansion-distance`: Either an integer or a list of integers. Each value is treated as an independent segmentation method, which is named as "10x_mm\_\_{value}um".

   - `localcores`: The maximum number of threads to use. Either an integer or a list of integers. In the former case, it will be broadcasted into a list if the value of `expansion-distance` is a list.

   - `localmem`: The maximum amount of memory (in GB) to use. Similar to `localcores`.

   - `_other_options`: For options other than those above, in a form similar to: "--option_1 value_1 --option_2 value_2". Note that options `--boundary-stain` and `--interior-stain` provided by the user here will be filtered out.

2. `10x`
   10X Xenium Ranger's cell segmentation algorithm is used, where both `--boundary-stain` and `--interior-stain` are disabled, so only DAPI signal is used.

   - `expansion-distance`: Similar to above. Each value is treated as an independent segmentation method, which is named as "10x\_{value}um".

   - `localcores`: Similar to above.

   - `localmem`: Similar to above.

   - `_other_options`: Similar to above.

3. `baysor`
   [Baysor (v0.7.0+)](https://github.com/kharchenkolab/Baysor) is used for segmentation. Users can specify some options of the method using the following keys:

   - `_config`: The path, either absolute or relative, to the configuration file in TOML format for Baysor. If a relative path is given, it must be relative to the current working directory. By default the recommended configuration (stored in `workflow/configs/baysor_xenium.toml`) by developers of Baysor is used, assuming the current working directory is the root of this repo.

   - `_threads`: The maximum number of threads to use.

   - `_memory`: The maximum amount of memory (in GB) to use.

   - `_other_options`: For options not defined in the configuration, in a form similar to that mentioned above.

4. `proseg`
   [Proseg (v2.0.0+)](https://github.com/dcjones/proseg) is used for segmentation. Users can specify some options of the method using the following keys:

   - `_threads`: The maximum number of threads to use.

   - `_memory`: The maximum amount of memory (in GB) to use.

   - `_other_options`: For options not defined in the configuration, in a form similar to that mentioned above.

5. `segger`
   [Segger](https://github.com/EliHei2/segger_dev) is used for segmentation. Users can specify some options at different stages during the computation.

   - `preprocess` stage

     - `tile_width`: Width of the tiles in pixels. Ignored if `tile_size` is provided. `200` by default.
     - `tile_height`: height of the tiles in pixels. Ignored if `tile_size` is provided. `200` by default.
     - `_threads`: Number of workers for parallel processing. `1` by default.
     - `_other_options`: For options not defined in the configuration except for `--sample_type`, in a form similar to that mentioned above.

   - `train` stage

     - `num_tx_tokens`: Number of tokens to encode genes in the panel. `500` by default. The value should not be smaller than the number of genes in the panel. For "5k" gene panel, the value is recommended to be set to `6000`.
     - `accelerator`: Device type to use for training (e.g.,`cuda`, `cpu`). `cpu` by default.
     - `devices`: Number of devices (GPUs) to use. `0` by default.
     - `max_epochs`: Number of epochs for training. `200` by default.
     - `batch_size`: Batch size for training. `4` by default.
     - `_threads`: Number of workers for data loading. `2` by default.
     - `_other_options`: For options not defined in the configuration except for `--sample_tag`, in a form similar to that mentioned above.

   - `predict` stage

     - `use_cc`: Use connected components if specified. `False` by default.
     - `_other_options`: For options not defined in the configuration except for `--batch_size`, `--knn_method`, and `--cell_id_col`, in a form similar to that mentioned above.

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

2. `normalisation`
   Normalisation methods to use during standard Seurat analysis, defined under `methods`:

   - `lognorm`: Log-normalisation.
   - `sctransform`: sctransform.

   At least of one of the two methods should be used when running the workflow.

3. `dim_reduction`
   Parameters used in dimension reduction.

   - `n_dims`: Number of dimensions to keep. `50` by default.

4. `clustering`
   Parameters used in clustering.

   - `resolution`: Parameter for `Seurat::FindClusters`. `0.8` by default.

## Coexpression (key: `coexpression`)

This section specifies a list of methods (with key `methods`) and a list of target counts (with key `target_counts`) for computing coexpression per segmentation method per sample.

1. `methods`
   The workflow supports the following methods: `conditional`, `jaccard`, `pearson`, and `spearman`.

2. `target_counts`
   By default the workflow uses `30`, `50`, and `200` globally for all samples. Users can specify different values for certain panels in `experiments.yml`.

## Cell type annotation (key: `cell_type_annotation`)

This section for now only supports reference based approach (specified with key `reference_based`) for cell type annotation.

1. `reference_based`
   Use reference-based methods for cell type annotation, where condition-specific single-cell RNA sequencing datasets are used as references (please refer to the details in section `experiments.yml` above).

   With key `methods`, this section specifies a number of methods that could be used for cell type annotation, including `rctd_class_aware`, `rctd_class_unaware`, `singler`, `seurat`, and `xgboost`. Command line arguments for each of the methods are configured in a separate dictionary with the method name being the key. Removing any methods inside `methods` will rule them out from the workflow.

   For both `rctd_class_aware` and `rctd_class_unaware`, RCTD is used for the annotation. The difference between them is that for the former the parent class of the current one in the reference is passed as an argument, while for the latter it is agnostic about the parent class.

   Users can use another key `modes` to specify the the mode under which to run annotation methods. For now only `single_cell` mode is supported.

   - `rctd`
     RCTD is used for annotation.

     - `UMI_min_sigma`: #TODO
     - `_threads`: The maximum amount of memory (in GB) to use. `20` bu default.
     - `_other_options`: For options not defined in the configuration, in a form similar to that mentioned above.

   - `singler`
     SingleR is used for annotation.

     - `genes`: #TODO
     - `de_method`: #TODO
     - `aggr_ref`: #TODO
     - `aggr_args_rank`: #TODO
     - `aggr_args_power`: #TODO
     - `_other_options`: Same as above.

   - `seurat`
     Seurat is used for annotation.

     - `min_dim`: #TODO
     - `max_dim`: #TODO
     - `_other_options`: Same as above.

   - `xgboost`
     XGboost is used for annotation.

     - `nrounds`: #TODO
     - `eta`: #TODO
     - `_other_options`: Same as above.
