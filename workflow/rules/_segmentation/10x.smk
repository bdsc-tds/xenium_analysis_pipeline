#######################################
#                Rules                #
#######################################

rule run10x:
    input:
        data_dir=f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{{sample_id}}'
    output:
        directory(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/raw_results')
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/run10x.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}',
        abs_input=lambda wildcards, input: os.path.abspath(
            input["data_dir"]
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/logs/run10x.log'
        ),
        expansion_distance=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.segmentation_id,
            "expansion-distance"
        ),
        localmem=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.segmentation_id,
            "localmem"
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.segmentation_id,
            "_other_options",
            replace_none=""
        )
    threads:
        lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.segmentation_id,
            "localcores"
        )
    resources:
        mem_mb=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.segmentation_id,
            "localmem"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger resegment --id=raw_results "
        "--xenium-bundle={params.abs_input} "
        "--expansion-distance={params.expansion_distance} "
        "--localcores={threads} "
        "--localmem={params.localmem} "
        "{params.other_options} &> {params.abs_log}"
