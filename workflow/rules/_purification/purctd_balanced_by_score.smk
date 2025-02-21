#######################################
#                Rules                #
#######################################

rule runPuRCTDBalancedByScore:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds', # annotation_id to be replaced with "matched_reference/rctd_class_aware/Leveli/single_cell"
        spatial_neighborhood_scores=f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/spatial_neighborhood_scores.parquet' ## annotation_id to be replaced with "matched_reference/rctd_class_aware/Leveli/single_cell"
        fully_purified_counts=protected(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/fully_purified_counts.h5'),
        fully_based_purified_counts_metadata=(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/fully_purified_counts_metadata.parquet')
    output:
        score_purified_counts=protected(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/score_purified_counts.h5'),
        score_purified_counts_metadata=protected(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/score_purified_counts_metadata.parquet')
    params:
        score_name="neighborhood_weights_second_type"
        score_threshold=0.1
    log:
        f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/logs/runPuRCTDBalancedBySpotClass.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_purification/purctd_balanced_by_score.R"
