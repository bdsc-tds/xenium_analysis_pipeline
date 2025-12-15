#######################################
#              Functions              #
#######################################

def get_generated_cells_geojson(wildcards) -> str:
    """
    Cells GeoJSON generated upstream by the pipeline.
    Convention:
      {output_path}/custom_segmentation/{compact_segmentation_id}/{sample_id}/cells.geojson
    """
    p = f'{config["output_path"]}/custom_segmentation/{wildcards.compact_segmentation_id}/{wildcards.sample_id}/cell_boundaries.geojson'
    return normalise_path(p, return_dir=False, check_exist=True)


#######################################
#                Rules                #
#######################################

rule generateVisiumGrid:
    input:
        xenium_bundle=get_input2_or_params4run10x,
    output:
        cells=get_generated_cells_geojson
    log:
        f'{config["output_path"]}/custom_segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/generateVisiumGrid.log'
    params:
      diameter=55, # ideally defined from the segmentation name... eg., grid_visium_55um 
      
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_custom_segmentation/_grid/visium_grid.R"
      
