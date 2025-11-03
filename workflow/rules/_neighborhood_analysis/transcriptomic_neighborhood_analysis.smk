#######################################
#                Rules                #
#######################################

rule runTranscriptomicNeighborhoodAnalysis:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        transcriptomic_neighborhood_scores=protected(f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/transcriptomic_neighborhood_scores.parquet')
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        DO_prune=False,
        k_knn=10
    log:
        f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runTranscriptomicNeighborhoodAnalysis.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        annotation_id=r".+/rctd_.+"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt * 50,
            51200,
        )
    script:
        "../../scripts/_neighborhood_analysis/transcriptomic_neighborhood_analysis.R"

use rule runTranscriptomicNeighborhoodAnalysis as runPostCountCorrectionBySplitTranscriptomicNeighborhoodAnalysis with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        transcriptomic_neighborhood_scores=protected(f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/transcriptomic_neighborhood_scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionBySplitTranscriptomicNeighborhoodAnalysis.log'
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r".+/rctd_.+"

use rule runTranscriptomicNeighborhoodAnalysis as runPostCountCorrectionByOvrlpyTranscriptomicNeighborhoodAnalysis with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        transcriptomic_neighborhood_scores=protected(f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/transcriptomic_neighborhood_scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByOvrlpyTranscriptomicNeighborhoodAnalysis.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r".+/rctd_.+"

use rule runTranscriptomicNeighborhoodAnalysis as runPostCountCorrectionByResolviUnsupervisedTranscriptomicNeighborhoodAnalysis with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        transcriptomic_neighborhood_scores=protected(f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/transcriptomic_neighborhood_scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviUnsupervisedTranscriptomicNeighborhoodAnalysis.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r".+/rctd_.+"

use rule runTranscriptomicNeighborhoodAnalysis as runPostCountCorrectionByResolviSupervisedTranscriptomicNeighborhoodAnalysis with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        transcriptomic_neighborhood_scores=protected(f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/transcriptomic_neighborhood_scores.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/logs/runPostCountCorrectionByResolviSupervisedTranscriptomicNeighborhoodAnalysis.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r".+/rctd_.+"
