#######################################
#                Rules                #
#######################################

rule runScdblfinder:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        doublet_scores=protected(f'{config["output_path"]}/doublet_finding/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{doublet_id}}/doublet_scores.parquet')
    params:
        label_aware=lambda wildcards: re.fullmatch(
            "scdblfinder_label_aware",
            wildcards.doublet_id,
        ) is not None
    log:
        f'{config["output_path"]}/doublet_finding/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{doublet_id}}/logs/runScdblfinder.log'
    wildcard_constraints:
        annotation_id=r".+/rctd_.+",
        doublet_id=r"scdblfinder_.+"
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt * 100,
            1024000,
        )
    script:
        "../../scripts/_doublet_finding/scdblfinder.R"
