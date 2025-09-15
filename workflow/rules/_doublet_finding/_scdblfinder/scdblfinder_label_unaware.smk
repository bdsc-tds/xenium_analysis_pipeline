#######################################
#                Rules                #
#######################################

rule runScdblfinderLabelUnaware:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        doublet_scores=protected(f'{config["output_path"]}/doublet_finding/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scdblfinder_label_unaware/doublet_scores.parquet')
    log:
        f'{config["output_path"]}/doublet_finding/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scdblfinder_label_unaware/logs/runScdblfinderLabelUnaware.log'
    wildcard_constraints:
        annotation_id=r".+/rctd_.+"
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt**2 * 100,
            1024000,
        )
    script:
        "../../../scripts/_doublet_finding/scdblfinder_label_unaware.R"
