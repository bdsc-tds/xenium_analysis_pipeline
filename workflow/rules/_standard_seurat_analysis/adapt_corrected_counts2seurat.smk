#######################################
#              Functions              #
#######################################




#######################################
#                Rules                #
#######################################

rule adaptCorrectedCounts2Seurat:
    input:
        raw_obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds',
        corrected_counts=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/corrected_counts.h5'
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/raw_seurat.rds')
    params:

    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/logs/adaptCorrectedCounts2Seurat.log'
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**2 * 2048, 512000)
    script:
        "../../scripts/_standard_seurat_analysis/adapt_corrected_counts2seurat.R"
