#######################################
#                Rules                #
#######################################

rule runStandardSCTdimredClust:
    input:
        xe=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/lognorm_seurat.rds'
    output:
        xe=protected(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/preprocessed_seurat.rds') 
    params:
        default_assay=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "default_xenium_assay"
        ),
        n_dims=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "dim-reduction",
            "n_dims"
        ),
        resolution=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "clustering",
            "resolution"
        )
    log:
        f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/out' # @Senbai, can/should we use the same logfile?
    script:
        '../../scripts/_standard_seurat_analysis/standard_QC_xenium.R'  
