#######################################
#                Rules                #
#######################################

rule runStandardSCTdimredClust:
    input:
        xe_raw=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/raw_seurat.rds'
        xe=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/preprocessed_seurat.rds'
    output:
        report=protected(f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/reports/standard_seurat_analysis.html')
    params:
        sample_id=f'{{sample_id}}', # Not sure
        default_assay=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "default_xenium_assay"
        )
    log:
        f'{config["outputDir"]}/segmentation/{{segmentation_id}}/{{sample_id}}/reports/out' # @Senbai, can/should wwe use the same logfile?
    script:
        '../../scripts/_standard_seurat_analysis/generate_report.R'  
