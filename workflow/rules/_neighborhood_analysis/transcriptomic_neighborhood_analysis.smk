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
        DO_prune=False, 
        k_knn=20
    log:
        f'{config["output_path"]}/neighborhood_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runTranscriptomicNeighborhoodAnalysis.log'
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(
            input.size_mb * attempt * 10,
            1024,
        )
    script:
        "../../scripts/_neighborhood_analysis/transcriptomic_neighborhood_analysis.R"
