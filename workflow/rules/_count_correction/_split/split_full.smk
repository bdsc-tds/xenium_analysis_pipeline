#######################################
#                Rules                #
#######################################

rule runSplitFull:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        corrected_counts=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/split_full/corrected_counts.h5'),
        corrected_counts_metadata=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/split_full/corrected_counts_metadata.parquet')
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/split_full/logs/runSplitFull.log'
    wildcard_constraints:
        annotation_id=r".+/rctd_.+"
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt**2 * 100, #70
            1024000,
        )
    script:
        "../../../scripts/_count_correction/split_full.R"
