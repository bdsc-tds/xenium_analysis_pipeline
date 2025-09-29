#######################################
#              Functions              #
#######################################

def get_rctd_results4runReferenceBasedTACIT(wildcards) -> str:
    prefix = f'{config["output_path"]}/cell_type_annotation/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}'
    rctd_id = re.sub(
        "(?<=/)tacit(?=/)",
        "rctd_class_aware",
        wildcards.annotation_id,
    )
    return f"{prefix}/{rctd_id}/output.rds"


#######################################
#                Rules                #
#######################################

rule runReferenceBasedTACIT:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation,
        rctd=get_rctd_results4runReferenceBasedTACIT
    output:
        tacit_output=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        tacit_weights=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
    params:
        mode=lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
            "tacit",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "mode",
        ),
        N_binary_markers_per_cell_type = 10,
        N_PCs = 50,
        r = 10, # minicluster size in TACIT
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
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/tacit/.+"
    threads:
        lambda wildcards: max(
            get_dict_value(
                config,
                "cell_type_annotation",
                extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
                "tacit",
                extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
                "_threads",
                replace_none=2,
            ),
            2,
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 20, 20480)
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedTACIT.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/tacit.R"

# use rule runReferenceBasedTACIT as runPostCountCorrectionBySplitReferenceBasedTACIT with:
#     input:
#         query=lambda wildcards: get_path2query4annotation(
#             wildcards,
#             is_post_correction=True,
#         ),
#         reference=get_path2reference4reference_based_annotation,
#         rctd=get_rctd_results4runReferenceBasedTACIT
#     output:
#         tacit_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
#         labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
#         tacit_weights=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
#     log:
#         f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionBySplitReferenceBasedTACIT.log'
#     wildcard_constraints:
#         count_correction_id=r"split.+",
#         annotation_id=r"reference_based/.+/tacit/.+"

use rule runReferenceBasedTACIT as runPostCountCorrectionByOvrlpyReferenceBasedTACIT with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation,
        rctd=get_rctd_results4runReferenceBasedTACIT
    output:
        tacit_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        tacit_weights=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByOvrlpyReferenceBasedTACIT.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/tacit/.+"

use rule runReferenceBasedTACIT as runPostCountCorrectionByResolviUnsupervisedReferenceBasedTACIT with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation,
        rctd=get_rctd_results4runReferenceBasedTACIT
    output:
        tacit_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        tacit_weights=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviUnsupervisedReferenceBasedTACIT.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/tacit/.+"

use rule runReferenceBasedTACIT as runPostCountCorrectionByResolviSupervisedReferenceBasedTACIT with:
    input:
        query=lambda wildcards: get_path2query4annotation(
            wildcards,
            is_post_correction=True,
        ),
        reference=get_path2reference4reference_based_annotation,
        rctd=get_rctd_results4runReferenceBasedTACIT
    output:
        tacit_output=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        tacit_weights=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviSupervisedReferenceBasedTACIT.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r"reference_based/.+/tacit/.+"

