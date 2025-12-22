#######################################
#              Functions              #
#######################################

def get_generated_cells_geojson(wildcards) -> str:
    """
    Cells GeoJSON generated upstream by the pipeline.
    Convention:
      {output_path}/custom_segmentation/{compact_segmentation_id}/{sample_id}/cells.geojson
    """
    p = f'{config["output_path"]}/segmentation/{wildcards.compact_segmentation_id}/{wildcards.sample_id}/processed_results/cell_boundaries.geojson'
    return p


def get_input2_or_params4generateVisiumGrid(wildcards, for_input: bool = True) -> str:
    use_raw_data, ret = get_raw_data_dir(wildcards.sample_id)

    if for_input:
        return ret

    if use_raw_data:
        return normalise_path(
            ret,
            pat_flags=re.IGNORECASE,
            return_dir=True,
            check_exist=True,
        )

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False,
    )


#######################################
#                Rules                #
#######################################

rule generateVisiumGrid:
    input:
        xenium_bundle=get_input2_or_params4generateVisiumGrid
    output:
        cells=f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/processed_results/cell_boundaries.geojson'
    log:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/generateVisiumGrid.log'
    params:
      diameter=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "diameter"
        ), # 55 by default and in most cases
    wildcard_constraints:
        compact_segmentation_id=r"grid_visium_\w+um" # should be grid_visium_55um
    resources: # not demanding reads cell_boundaries.geojson from the xenium bundles and create another one. 
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_segmentation/_grid/visium_grid.R"
      
