#######################################
#                Rules                #
#######################################

rule runSpatialNeighborhoodAnalysis:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        post_processed_rctd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'
    output:
        spatial_neighborhood_scores=protected(f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/spatial_neighborhood_scores.parquet')
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        DO_prune=True,
        pruning_radius=15,
        k_knn=20
    log:
        f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runSpatialNeighborhoodAnalysis.log'
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
        "../../scripts/_neighborhood_analysis/spatial_neighborhood_analysis.R"
