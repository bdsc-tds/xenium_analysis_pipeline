#######################################
#                Rules                #
#######################################

rule runReferenceBasedSeurat:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet')
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        annotation_id=lambda wildcards: wildcards.annotation_id,
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
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 10240),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedSeurat.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/seurat.R"

use rule runReferenceBasedSeurat as runPostCountCorrectionReferenceBasedSeurat with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionReferenceBasedRCTD.log'
    wildcard_constraints:
        count_correction_id=r"(?!ovrlpy)(?!resolvi_unsupervised)",
        annotation_id=r"reference_based/.+/seurat/.+"

use rule runReferenceBasedSeurat as runPostCountCorrectionByOvrlpyReferenceBasedSeurat with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByOvrlpyReferenceBasedSeurat.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/seurat/.+"

use rule runReferenceBasedSeurat as runPostCountCorrectionByResolviUnsupervisedReferenceBasedSeurat with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation
    output:
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviUnsupervisedReferenceBasedSeurat.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/seurat/.+"
