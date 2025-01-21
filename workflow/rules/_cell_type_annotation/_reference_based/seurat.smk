#######################################
#                Rules                #
#######################################

rule runReferenceBasedSeurat:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/scores.parquet')
    params:
        future_globals_maxSize=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "_future_globals_maxSize",
            replace_none=1,
        ) * 10**9,
        annotation_id=lambda wildcards: wildcards.annotation_id,
        ref_default_assay=cac.REF_SEURAT_DEFAULT_ASSAY,
        xe_default_assay=cac.XE_SEURAT_DEFAULT_ASSAY,
        REF_MIN_UMI=cac.REF_MIN_UMI,
        REF_MAX_UMI=cac.REF_MAX_UMI,
        XE_MIN_UMI=cac.XE_MIN_UMI,
        XE_MIN_counts=cac.XE_MIN_counts,
        cell_min_instance=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_CELL_MIN_INSTANCES_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0])[0],
            extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
            replace_none=25,
        ),
        min_dim=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "seurat",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "min_dim",
            replace_none=1,
        ),
        max_dim=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "seurat",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "max_dim",
            replace_none=50,
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/seurat/.+"
    resources:
        mem_mb=lambda wildcards, input: max(input.size_mb * 30, 10240)
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/logs/runReferenceBasedSeurat.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/seurat.R"
