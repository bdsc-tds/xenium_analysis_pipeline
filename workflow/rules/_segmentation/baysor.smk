#######################################
#                Rules                #
#######################################

rule runBaysor:
    input:
        data=f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{{sample_id}}/transcripts.parquet'
    output:
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_2d.json'),
        protected(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_3d.json')
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/runBaysor.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results',
        abs_input=lambda wildcards, input: os.path.abspath(
            input["data"]
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
        raw_data=f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{{sample_id}}',
        segmentation=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation.csv',
        polygons=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/raw_results/segmentation_polygons_3d.json' #TODO: 2d or 3d?
    output:
        directory(f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/baysor/{{sample_id}}/logs/normaliseBaysor.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/baysor/{{sample_id}}',
        abs_input_raw_data=lambda wildcards, input: os.path.abspath(
            input["raw_data"]
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
            "localmem"
        )
    threads:
        get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "localcores"
        )
    resources:
        mem_mb=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "localmem"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger import-segmentation --id=normalised_results "
        "--xenium-bundle {params.abs_input_raw_data} "
        "--transcript-assignment={params.abs_input_segmentation} "
        "--viz-polygons={params.abs_input_polygons} "
        "--units=microns "
        "--localcores={threads} "
        "--localmem={params.localmem} &> {params.abs_log}"
