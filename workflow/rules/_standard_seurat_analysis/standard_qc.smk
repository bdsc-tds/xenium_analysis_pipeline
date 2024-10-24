#######################################
#                Rules                #
#######################################

rule runStandardQC:
    input:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/raw_seurat.rds'
    output:
        temp(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/qced_seurat.rds')
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
                "min_counts",
                replace_none=10
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
                "min_features",
                replace_none=5
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
                "max_counts",
                replace_none=float("inf")
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
                "max_features",
                replace_none=float("inf")
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
                "min_cells",
                replace_none=1
            )
        ),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        default_layer=sec.SEURAT_DEFAULT_LAYER
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/runStandardQC.log'
    script:
        "workflow/scripts/_standard_seurat_analysis/standard_qc.R"
