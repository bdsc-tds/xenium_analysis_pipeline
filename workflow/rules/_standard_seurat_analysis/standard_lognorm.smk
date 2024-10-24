#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/qced_seurat.rds'
    output:
        temp(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/lognormed_seurat.rds')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/runStandardLogNorm.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_lognorm.R"
