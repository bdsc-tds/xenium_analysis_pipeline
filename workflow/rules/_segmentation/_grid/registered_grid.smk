#######################################
#              Functions              #
#######################################

def get_registered_cells_geojson(wildcards) -> str:
    """
    Cells or grid units GeoJSON [manually] registered upstream.
    """
    p = f'{config["output_path"]}/registration/{wildcards.compact_segmentation_id}/{wildcards.sample_id}/cell_boundaries.geojson'
    return p

#######################################
#                Rules                #
#######################################

rule importRegisteredGrid:
    input:
        registered_cells=get_registered_cells_geojson
    output:
        cells=f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/processed_results/cell_boundaries.geojson'
    log:
        f'{config["output_path"]}/registration/{{compact_segmentation_id}}/{{sample_id}}/logs/importRegisteredGrid.log'
    wildcard_constraints:
        compact_segmentation_id=r"grid_registered_\w+um" # should be grid_registered_visium_55um
    resources: # not demanding reads cell_boundaries.geojson from the xenium bundles and create another one. 
    shell:
        r"""
        mkdir -p $(dirname {output.cells})
        ln -sf $(realpath {input.registered_cells}) {output.cells}
        echo "Linked {input.registered_cells} to {output.cells}" > {log}
        """
      
