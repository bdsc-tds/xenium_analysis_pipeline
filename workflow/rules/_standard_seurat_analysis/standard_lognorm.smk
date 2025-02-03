#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_seurat.rds'
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        normalisation_id="lognorm"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 50, 10240)
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/logs/runStandardLogNorm.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_lognorm.R"
