#######################################
#              Functions              #
#######################################

def get_input2_or_params4run10x(wildcards, for_input: bool = True) -> str:
    gene_panel_file: str | None = get_gene_panel_file(wildcards.sample_id, config)

    with open(checkpoints.check10xVersions.get(sample_id=wildcards.sample_id).output[0], "r", encoding="utf-8") as fh:
        versions: dict[str, Any] = json.load(fh)
    
    matched: bool = get_dict_value(
        versions,
        "match",
        str(get_dict_value(
            config,
            "reprocess",
            "level"
        ))
    )

    use_raw_data: bool = True if gene_panel_file is None or matched else False
    if use_raw_data:
        ret: str = f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{wildcards.sample_id}'
    else:
        ret = f'{config["output_path"]}/reprocessed/{wildcards.sample_id}/results'

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
        directory(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/run10x.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}',
        abs_input=lambda wildcards: os.path.abspath(
            get_input2_or_params4run10x(wildcards, for_input=False)
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
    wildcard_constraints:
        segmentation_id=r"10x_\w*?_?\d+?um"
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
        "xeniumranger resegment --id=normalised_results "
        "--xenium-bundle={params.abs_input} "
        "--expansion-distance={params.expansion_distance} "
        "--localcores={threads} "
        "--localmem={params.localmem} "
        "{params.other_options} &> {params.abs_log}"
