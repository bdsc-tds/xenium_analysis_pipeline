#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        xe=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/QC_seurat.rds'
    output:
        xe=temp(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/lognorm_seurat.rds') # @Senbai, should I use the same name for all the temp objects or its better to use "qc_seurat", "longnorm_seurat", "sctransform_seurat" etc?
    params:
      default_assay=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "default_xenium_assay"
        )
    log:
        f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/out' # @Senbai, can/should we use the same logfile?
    script:
        '../../scripts/_standard_seurat_analysis/standard_lognorm_xenium.R'  
