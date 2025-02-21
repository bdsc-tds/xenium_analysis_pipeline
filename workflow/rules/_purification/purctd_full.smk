#######################################
#                Rules                #
#######################################

rule runPuRCTDFullPurification:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds' # annotation_id to be replaced with "matched_reference/rctd_class_aware/Leveli/single_cell"
    output:
        fully_purified_counts=protected(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/fully_purified_counts.h5'),
        fully_based_purified_counts_metadata=protected(f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/fully_purified_counts_metadata.parquet')
    params:
    log:
        f'{config["output_path"]}/purification/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/purctd/{{reference_name}}/logs/runPuRCTDFullPurification.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_purification/purctd_full.R"
