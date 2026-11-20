#######################################
#              Functions              #
#######################################



#######################################
#                Rules                #
#######################################

rule generateSquareGrid:
    input:
        xenium_bundle=get_input2_or_params4generateVisiumGrid
    output:
        cells=f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/processed_results/cell_boundaries.geojson'
    log:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/generateSquareGrid.log'
    params:
      bin_size=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "bin_size"
        ), # 2, 8, 16 by default and in most cases
    wildcard_constraints:
        compact_segmentation_id=r"grid_square_\w+um" # should be grid_square_2um or better grid_square_002um
    resources: # not demanding reads cell_boundaries.geojson from the xenium bundles and create another one. 
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_segmentation/_grid/square_grid.R"
      
