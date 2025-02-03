#######################################
#                Rules                #
#######################################

rule runReferenceBasedXGBoost:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet')
    params:
        annotation_id=lambda wildcards: wildcards.annotation_id,
        ref_assay=lambda wildcards: get_assay_name(wildcards, True),
        xe_assay=lambda wildcards: get_assay_name(wildcards, False),
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
        nrounds=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "xgboost",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "nrounds",
            replace_none=1000,
        ),
        eta=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "xgboost",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "eta",
            replace_none=0.3,
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/xgboost/.+"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 50, 10240)
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedXGBoost.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/xgboost.R"
