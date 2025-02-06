#######################################
#              Functions              #
#######################################

def get_input2_or_params4run10x(wildcards, for_input: bool = True) -> str:
    use_raw_data, ret = get_raw_data_dir(wildcards.sample_id)

    if for_input:
        return ret

    if use_raw_data:
        return normalise_path(
            ret,
            pat_flags=re.IGNORECASE,
            return_dir=True,
            check_exist=True
        )

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False
    )


#######################################
#                Rules                #
#######################################

rule run10x:
    input:
        get_input2_or_params4run10x
    output:
        directory(f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/run10x.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}',
        abs_input=lambda wildcards: os.path.abspath(
            get_input2_or_params4run10x(wildcards, for_input=False)
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/{wildcards.compact_segmentation_id}/{wildcards.sample_id}/logs/run10x.log'
        ),
        expansion_distance=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "expansion-distance"
        ),
        localmem=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "localmem"
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "_other_options",
            replace_none=""
        )
    wildcard_constraints:
        compact_segmentation_id=r"10x_\w*?_?\d+?um"
    threads:
        lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "localcores"
        )
    retries:
        0
    resources:
        mem_mb=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            wildcards.compact_segmentation_id,
            "localmem"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger resegment --id=normalised_results "
        "--xenium-bundle={params.abs_input} "
        "--expansion-distance={params.expansion_distance} "
        "--localcores={threads} "
        "--localmem={params.localmem} "
        "{params.other_options} &> {params.abs_log}"
