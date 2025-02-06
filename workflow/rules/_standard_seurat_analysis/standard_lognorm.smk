#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=temp(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/normalised_seurat.rds'),
        cells=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/cells.parquet'),
        counts=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/counts.parquet')
    params:
        normalised_assay=sec.SEURAT_DEFAULT_ASSAY,
        normalised_layer=sec.SEURAT_ALT_LAYER,
        normalisation_id="lognorm"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt**2 * 10, 10240)
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/logs/runStandardLogNorm.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_lognorm.R"
