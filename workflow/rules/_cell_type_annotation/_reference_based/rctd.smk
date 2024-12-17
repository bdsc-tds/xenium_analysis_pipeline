#######################################
#              Functions              #
#######################################

def get_class_level4runReferenceBasedRCTD(wildcards) -> str:
    current_level: str = extract_layers_from_experiments(wildcards.annotation_id, [3])[0]
    available_levels: list[str] = get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_LEVELS_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0])[0],
            extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
        )

    try:
        current_level_idx: int = available_levels.index(current_level)
    except ValueError:
        raise RuntimeError(f"Error! Cannot find {current_level} in [{",".join(available_levels)}].")
    
    return available_levels[0] if current_level_idx == 0 else available_levels[current_level_idx - 1]


#######################################
#                Rules                #
#######################################

rule runReferenceBasedRCTD:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/labels.csv'),
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/scores.csv')
    params:
        annotation_id=lambda wildcards: wildcards.annotation_id,
        class_level=get_class_level4runReferenceBasedRCTD,
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
        UMI_min_sigma=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "UMI_min_sigma",
        ),
        cores=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            "reference_based",
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "_threads",
            replace_none=20,
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd/.+"
    threads:
        lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            "reference_based",
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "_threads",
            replace_none=20,
        )
    resources:
        mem_mb=lambda wildcards, input: max(input.size_mb * 30, 20480)
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/logs/runReferenceBasedRCTD.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/rctd.R"
