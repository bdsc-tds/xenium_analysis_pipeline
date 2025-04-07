#######################################
#              Functions              #
#######################################

def get_class_level4runReferenceBasedRCTD(wildcards) -> str | None:
    method_name: str = extract_layers_from_experiments(wildcards.annotation_id, [2])[0]

    if method_name == "rctd_class_unaware":
        return None
    elif method_name == "rctd_class_aware":
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
    else:
        raise RuntimeError(f"Error! Unknown method name: {method_name}.")


def is_int_counts4runReferenceBasedRCTD(wildcards) -> bool:
    return False if re.match(
        "proseg_expected",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) else True


#######################################
#                Rules                #
#######################################

rule runReferenceBasedRCTD:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        rds_output=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        scores=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet'),
        out_res_df=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet'),
        out_w=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights.parquet'),
        out_wd=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights_doublet.parquet'),
        out_sc=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/singlet_score.parquet'),
        out_sm=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/score_mat.parquet')
    params:
        annotation_id=lambda wildcards: wildcards.annotation_id,
        class_level=get_class_level4runReferenceBasedRCTD,
        ref_assay=lambda wildcards: get_assay_name4annotation(wildcards, True),
        xe_assay=lambda wildcards: get_assay_name4annotation(wildcards, False),
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
        is_int_counts=is_int_counts4runReferenceBasedRCTD,
        cores=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "_threads",
            replace_none=20,
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd_.+"
    threads:
        lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "_threads",
            replace_none=20,
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 20, 20480)
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedRCTD.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/rctd.R"

use rule runReferenceBasedRCTD as runPostCountCorrectionReferenceBasedRCTD with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        rds_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        scores=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet'),
        out_res_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet'),
        out_w=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights.parquet'),
        out_wd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights_doublet.parquet'),
        out_sc=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/singlet_score.parquet'),
        out_sm=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/score_mat.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionReferenceBasedRCTD.log'
    wildcard_constraints:
        count_correction_id=COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT,
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTD as runPostCountCorrectionByOvrlpyReferenceBasedRCTD with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        rds_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        scores=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/scores.parquet'),
        out_res_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet'),
        out_w=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/weights.parquet'),
        out_wd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/weights_doublet.parquet'),
        out_sc=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/singlet_score.parquet'),
        out_sm=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/score_mat.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByOvrlpyReferenceBasedRCTD.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTD as runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTD with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        rds_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        scores=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet'),
        out_res_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet'),
        out_w=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights.parquet'),
        out_wd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights_doublet.parquet'),
        out_sc=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/singlet_score.parquet'),
        out_sm=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/score_mat.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTD.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/rctd_.+"
