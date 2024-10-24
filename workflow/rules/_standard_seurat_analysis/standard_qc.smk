#######################################
#                Rules                #
#######################################

rule runStandardQC:
    input:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/raw_seurat.rds'
    output:
        temp(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/qced_seurat.rds') # @Senbai, should I use the same name for all the temp objects or its better to use "qc_seurat", "longnorm_seurat", "sctransform_seurat" etc?
    params:
        min_counts=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_counts",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_counts"
            )
        ),
        min_features=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_features",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_features"
            )
        ),
        max_counts=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "max_counts",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "max_counts"
            )
        ),
        max_features=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "max_features",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "max_features"
            )
        ),
        min_cells=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_cells",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_cells"
            )
        ),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        default_layer=sec.SEURAT_DEFAULT_LAYER
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/runStandardQC.log'
    script:
        "workflow/scripts/_standard_seurat_analysis/standard_qc.R"
