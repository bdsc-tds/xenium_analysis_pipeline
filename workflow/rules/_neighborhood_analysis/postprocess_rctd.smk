#######################################
#                Rules                #
#######################################

rule runPostprocessRCTD:
    input:
        rctd_result=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostprocessRCTD.log'
    wildcard_constraints:
        annotation_id=r".+/rctd_.+"
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt * 50,
            51200,
        )
    script:
        "../../scripts/_neighborhood_analysis/postprocess_rctd.R"

use rule runPostprocessRCTD as runPostCountCorrectionBySplitPostprocessRCTD with:
    input:
        rctd_result=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionBySplitPostprocessRCTD.log'
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r".+/rctd_.+"

use rule runPostprocessRCTD as runPostCountCorrectionByOvrlpyPostprocessRCTD with:
    input:
        rctd_result=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByOvrlpyPostprocessRCTD.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r".+/rctd_.+"

use rule runPostprocessRCTD as runPostCountCorrectionByResolviUnsupervisedPostprocessRCTD with:
    input:
        rctd_result=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviUnsupervisedPostprocessRCTD.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r".+/rctd_.+"

use rule runPostprocessRCTD as runPostCountCorrectionByResolviSupervisedPostprocessRCTD with:
    input:
        rctd_result=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviSupervisedPostprocessRCTD.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r".+/rctd_.+"
