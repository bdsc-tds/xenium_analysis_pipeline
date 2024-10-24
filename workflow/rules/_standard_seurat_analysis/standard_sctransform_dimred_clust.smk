#######################################
#                Rules                #
#######################################

rule runStandardSctransformDimRedClust:
    input:
        # f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/lognorm_seurat.rds'
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/qced_seurat.rds'
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/preprocessed_seurat.rds') 
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        n_dims=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "dim_reduction",
            "n_dims",
            replace_none=50
        ),
        resolution=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "clustering",
            "resolution",
            replace_none=0.8
        )
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/runStandardSctransformDimRedClust.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_sctransform_dimred_clust.R"
