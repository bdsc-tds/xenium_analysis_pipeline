#######################################
#                Rules                #
#######################################

rule runStandardQC:
    input:
        xe=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/raw_seurat.rds'
    output:
        xe=temp(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/QC_seurat.rds') # @Senbai, should I use the same name for all the temp objects or its better to use "qc_seurat", "longnorm_seurat", "sctransform_seurat" etc?
    params:
        min_counts=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "QC",
            "min_counts"
        ),
        min_features=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "QC",
            "min_features"
        ),
        max_counts=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "QC",
            "min_counts"
        ),
        max_features=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "QC",
            "min_features"
        ),
        min_cells=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "QC",
            "min_cells"
        ),
        default_assay=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "default_xenium_assay"
        ),
        default_layer=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "default_xenium_layer"
        )
    log:
        f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/out' # @Senbai, can/should wwe use the same logfile?
    script:
        '../../scripts/_standard_seurat_analysis/standard_QC_xenium.R'  
