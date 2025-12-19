# This rule takes a grid as a geojson file and runs xeniumranger input-segmentation to generate 10x-normalized output for pseudo-visium
# The resuls are stored as an alternative segmentation under "segmentation" level 
# expected grids are: 
# - grid_visium_55um - for visium-like grid with 55um diameter spots
# - grid_registered_visium_55um - for the registered spots (visium registered to xenium adjuscent slice)
# - grid_binned_8um - for visiumHD-like 
# - grid_registered_binned_8um - for registered visiumHD if any 

# the grids like grid_visium_55um and grid_binned_8um are generated with the `_segmentation/_grid/` rules, 
# while the registered one are provided manually, but all of them are stored in `/segmentation/{grid_id == segmentation_id}/{sample_id}/cell_boundaries.geojson`

#######################################
#              Functions              #
#######################################



#######################################
#                Rules                #
#######################################

rule importGrid:
    input:
        xenium_bundle=get_input2_or_params4run10x, # can be raw segmentation from 10x  
        cells=get_generated_cells_geojson
    output:
        directory(f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/importGrid.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}',
        abs_bundle=lambda wc: os.path.abspath(get_input2_or_params4run10x(wc, for_input=False)),
        abs_cells=lambda wc: os.path.abspath(get_generated_cells_geojson(wc)),
        units=lambda wc: get_dict_value(config, "segmentation", wc.compact_segmentation_id, "units", replace_none="microns"),
        localmem=lambda wc: get_dict_value(config, "segmentation", wc.compact_segmentation_id, "localmem"),
        other_options=lambda wc: get_dict_value(config, "segmentation", wc.compact_segmentation_id, "_other_options", replace_none=""),
        abs_log=lambda wc: os.path.abspath(
            f'{config["output_path"]}/segmentation/{wc.compact_segmentation_id}/{wc.sample_id}/logs/import-segmentation.log'
        )
    wildcard_constraints:
        compact_segmentation_id=r"grid_(?:[a-zA-Z]+_)+\d+um" # can be grid_visium_55um or grid_registered_visium_55um or grid_binned_8um
    threads:
        lambda wc: get_dict_value(config, "segmentation", wc.compact_segmentation_id, "localcores")
    resources:
        mem_mb=lambda wc: get_dict_value(config, "segmentation", wc.compact_segmentation_id, "localmem") * 1024
    container:
        config["containers"]["10x"]
    shell:
        r"""
        mkdir -p {params.work_dir}/logs
        cd {params.work_dir}
        xeniumranger import-segmentation --id=normalised_results \
          --xenium-bundle={params.abs_bundle} \
          --cells={params.abs_cells} \
          --units={params.units} \
          --localcores={threads} \
          --localmem={params.localmem} \
          {params.other_options} &> {params.abs_log}
        """
