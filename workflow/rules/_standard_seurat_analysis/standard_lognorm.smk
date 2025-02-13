#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/scale_data.parquet')
    params:
        assay=sec.SEURAT_DEFAULT_ASSAY,
        data_layer=sec.SEURAT_DATA_LAYER,
        scale_data_layer=sec.SEURAT_SCALE_DATA_LAYER,
        normalisation_id="lognorm"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt**3 * 20, 10240)
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/logs/runStandardLogNorm.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_lognorm.R"
