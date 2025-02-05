#######################################
#                Rules                #
#######################################

rule runStandardDimRedClust:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/normalised_counts/normalised_seurat.rds'
    output:
        obj=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'),
        cells=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/cells.parquet'),
        pca=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/pca.parquet'),
        umap=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/umap.parquet')
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        n_dims=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "dim_reduction",
            "n_dims",
            replace_none=50,
        ),
        resolution=lambda wildcards: get_dict_value(
            config,
            "standard_seurat_analysis",
            "clustering",
            "resolution",
            replace_none=0.8,
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 20480),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/logs/runStandardDimRedClust.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_dimred_clust.R"
