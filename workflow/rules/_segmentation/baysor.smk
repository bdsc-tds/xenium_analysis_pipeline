#######################################
#              Functions              #
#######################################

def get_other_options4runBaysor(wildcards, input) -> str:
    ret: str = get_dict_value(
        config,
        "segmentation",
        "baysor",
        "_other_options",
        replace_none=""
    )

    with open(input["xr_version"], "r", encoding="utf-8") as fh:
        xr_version: int = get_dict_value(
            json.load(fh),
            "system_software_version",
            "0"
        )

    if xr_version < 4:
        ret = " ".join([ret, "--polygon-format=GeometryCollection"])

    return ret

def get_input2_or_params4normaliseBaysor(wildcards) -> dict[str, str]:
    ret: dict[str, str] = {
        "data_dir": get_input2_or_params4run10x(wildcards)
    }

    with open(checkpoints.check10xVersions.get(sample_id=wildcards.sample_id).output[0], "r", encoding="utf-8") as fh:
        xr_version: int = get_dict_value(
            json.load(fh),
            "system_software_version",
            "0"
        )
    
    if xr_version < 4:
        ret["segmentation"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/processed_results/segmentation.csv'
        ret["polygons"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/processed_results/segmentation_polygons_2d.json'
    else:
        ret["segmentation"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/raw_results/segmentation.csv'
        ret["polygons"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/raw_results/segmentation_polygons_2d.json'
    
    return ret


#######################################
#                Rules                #
#######################################

rule runBaysor:
    input:
        data_file=f'{config["output_path"]}/reprocessed/{{sample_id}}/transcripts_snappy.parquet',
        xr_version=f'{config["output_path"]}/reprocessed/{{sample_id}}/versions.json'
    output:
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_2d.json'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_3d.json')
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/runBaysor.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results',
        abs_input=lambda wildcards, input: os.path.abspath(input["data_file"]),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/logs/runBaysor.log'
        ),
        abs_config=lambda wildcards: os.path.abspath(
            get_dict_value(
                config,
                "segmentation",
                "baysor",
                "_config"
            )
        ),
        other_options=lambda wildcards, input: get_other_options4runBaysor(wildcards, input)
    threads:
        get_dict_value(
            config,
            "segmentation",
            "baysor",
            "_threads"
        )
    resources:
        mem_mb=get_dict_value(
            config,
            "segmentation",
            "baysor",
            "_memory"
        ) * 1024
    container:
        config["containers"]["baysor"]
    shell:
        "cd {params.work_dir} && "
        "JULIA_NUM_THREADS={threads} && "
        "baysor run -c {params.abs_config} "
        "{params.other_options} "
        "{params.abs_input} :cell_id &> {params.abs_log}"

rule adjustBaysorResults:
    input:
        segmentation=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv',
        polygons=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_2d.json'
    output:
        segmentation=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/processed_results/segmentation.csv',
        polygons=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/processed_results/segmentation_polygons_2d.json'
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/adjustBaysorResults.log'
    shell:
        "python3 workflow/scripts/adjust_baysor_results.py "
        "--inseg {input.segmentation} "
        "--inpoly {input.polygons} "
        "-l {log} "
        "--outseg {output.segmentation} "
        "--outpoly {output.polygons}"

rule normaliseBaysor:
    input:
        unpack(get_input2_or_params4normaliseBaysor)
    output:
        directory(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/normaliseBaysor.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}',
        abs_input_data_dir=lambda wildcards: os.path.abspath(
            get_input2_or_params4run10x(wildcards, for_input=False)
        ),
        abs_input_segmentation=lambda wildcards, input: os.path.abspath(
            get_input2_or_params4normaliseBaysor(wildcards)["segmentation"]
        ),
        abs_input_polygons=lambda wildcards, input: os.path.abspath(
            get_input2_or_params4normaliseBaysor(wildcards)["polygons"]
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/logs/normaliseBaysor.log'
        ),
        localmem=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_memory"
        )
    threads:
        get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_threads"
        )
    resources:
        mem_mb=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_memory"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger import-segmentation --id=normalised_results "
        "--xenium-bundle {params.abs_input_data_dir} "
        "--transcript-assignment={params.abs_input_segmentation} "
        "--viz-polygons={params.abs_input_polygons} "
        "--units=microns "
        "--localcores={threads} "
        "--localmem={params.localmem} &> {params.abs_log}"
