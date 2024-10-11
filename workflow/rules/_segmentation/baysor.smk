#######################################
#              Functions              #
#######################################

def get_input2_or_params4runBaysor(wildcards, for_input: bool = True) -> str:
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
            pat_anchor_file="transcripts.parquet",
            return_dir=False,
            check_exist=True
        )

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_anchor_file="transcripts.parquet",
        return_dir=False,
        check_exist=False
    )


#######################################
#                Rules                #
#######################################

rule runBaysor:
    input:
        get_input2_or_params4runBaysor
    output:
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_2d.json'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_3d.json')
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/runBaysor.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results',
        abs_input=lambda wildcards: os.path.abspath(
            get_input2_or_params4runBaysor(wildcards, for_input=False)
        ),
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
        other_options=get_dict_value(
            config,
            "segmentation",
            "baysor",
            "_other_options",
            replace_none=""
        )
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

rule normaliseBaysor:
    input:
        data_dir=lambda wildcards: get_input2_or_params4run10x(wildcards),
        segmentation=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv',
        polygons=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_3d.json' #TODO: 2d or 3d?
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
            input["segmentation"]
        ),
        abs_input_polygons=lambda wildcards, input: os.path.abspath(
            input["polygons"]
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
