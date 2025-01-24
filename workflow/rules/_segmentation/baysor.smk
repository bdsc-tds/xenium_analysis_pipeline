#######################################
#              Functions              #
#######################################

def _compare_version(ver: dict[str, int], against: tuple, is_lower_bound: bool) -> bool:
    assert len(against) > 0

    for idx, val in enumerate(against):
        tmp: int = get_dict_value(ver, str(idx))

        if tmp == val:
            continue
        
        if is_lower_bound:
            return tmp > val
        else:
            return tmp < val

    return True

def get_xeniumranger_version(
    version_file_path: str,
    min_version: tuple | None = None,
    max_version: tuple | None = None
) -> tuple | bool:
    assert min_version is None or isinstance(min_version, tuple)
    assert max_version is None or isinstance(max_version, tuple)

    with open(version_file_path, "r", encoding="utf-8") as fh:
        ver: dict[str, int] = get_dict_value(
            json.load(fh),
            "system_software_version"
        )

    meets_min, meets_max = None, None

    if min_version is not None and len(min_version) > 0:
        meets_min = _compare_version(ver, min_version, True)

    if max_version is not None and len(max_version) > 0:
        meets_max = _compare_version(ver, max_version, False)

    if meets_min is not None and meets_max is not None:
        return meets_min and meets_max
    elif meets_min is not None:
        return meets_min
    elif meets_max is not None:
        return meets_max
    else:
        return (v for _, v in ver.items())

def get_other_options4runBaysor(wildcards, input) -> str:
    ret: str = get_dict_value(
        config,
        "segmentation",
        "baysor",
        "_other_options",
        replace_none=""
    )

    meets_min: bool = get_xeniumranger_version(
        input["xr_version"],
        min_version=(10,)
    )

    if not meets_min:
        ret = " ".join([ret, "--polygon-format=GeometryCollection"])

    return ret

def get_input2_or_params4normaliseBaysor(wildcards) -> dict[str, str]:
    ret: dict[str, str] = {
        "data_dir": get_input2_or_params4run10x(wildcards)
    }

    meets_min: bool = get_xeniumranger_version(
        checkpoints.check10xVersions.get(sample_id=wildcards.sample_id).output[0],
        min_version=(3, 1)
    )

    if meets_min:
        ret["segmentation"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/raw_results/segmentation.csv'
        ret["polygons"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/raw_results/segmentation_polygons_2d.json'
    else:
        ret["segmentation"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/processed_results/segmentation.csv'
        ret["polygons"] = f'{config["output_path"]}/segmentation/baysor/{wildcards.sample_id}/processed_results/segmentation_polygons_2d.json'

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
    retries:
        RETRIES_NUM
    threads:
        get_dict_value(
            config,
            "segmentation",
            "baysor",
            "_threads"
        )
    resources:
        mem_mb=lambda wildcards, attempt: get_dict_value(
            config,
            "segmentation",
            "baysor",
            "_memory"
        ) * 1024 * attempt
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
    retries:
        RETRIES_NUM
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 2048)
    shell:
        "python3 workflow/scripts/_segmentation/adjust_baysor_results.py "
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
        abs_input_segmentation=lambda wildcards: os.path.abspath(
            get_input2_or_params4normaliseBaysor(wildcards)["segmentation"]
        ),
        abs_input_polygons=lambda wildcards: os.path.abspath(
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
