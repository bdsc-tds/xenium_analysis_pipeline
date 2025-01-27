#######################################
#                Rules                #
#######################################

rule runStandardSctransformDimRedClust:
    input:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/lognormed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/preprocessed_seurat.rds') 
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        n_dims=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "dim_reduction",
            "n_dims",
            replace_none=50,
        ),
        resolution=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "clustering",
            "resolution",
            replace_none=0.8,
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 100, 20480),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/runStandardSctransformDimRedClust.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_sctransform_dimred_clust.R"
