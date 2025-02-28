#######################################
#                Rules                #
#######################################

rule runPuRCTDScoreBalanced:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds',
        fully_purified_counts=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/purctd_full/corrected_counts.h5',
        fully_purified_counts_metadata=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/purctd_full/corrected_counts_metadata.parquet',
        spatial_neighborhood_scores=f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/spatial_neighborhood_scores.parquet'
    output:
        corrected_counts=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/purctd_score_balanced/corrected_counts.h5'),
        corrected_counts_metadata=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/purctd_score_balanced/corrected_counts_metadata.parquet')
    params:
        score_name="neighborhood_weights_second_type",
        score_threshold=0.1
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/purctd_score_balanced/logs/runPuRCTDScoreBalanced.log'
    wildcard_constraints:
        annotation_id=r".+/rctd_.+"
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(
            input.size_mb * attempt * 40,
            20480,
        )
    script:
        "../../../scripts/_count_correction/purctd_score_balanced.R"
