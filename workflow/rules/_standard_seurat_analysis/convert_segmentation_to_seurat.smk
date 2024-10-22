#######################################
#                Rules                #
#######################################

rule runSegmentationToSeurat:
    input:
        xr_dir=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/normalised_results'
    output:
        xe=protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/raw_seurat.rds')
    params:
      spatial_dimname=lambda wildcards: get_dict_value(
            config,
            "standard-seurat-analysis",
            "object",
            "spatial_dimname"
        )
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/objects/out'
    container: 
        config["containers"]["r"]
    script:
        '../../scripts/convert_segmentation_to_seurat.R' # Not sure about relative path 

