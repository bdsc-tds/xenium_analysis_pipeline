#######################################
#              Functions              #
#######################################

def get_input2_or_params4runProseg(wildcards, for_input: bool = True) -> str:
    use_raw_data, ret = get_raw_data_dir(wildcards.sample_id)

    if for_input:
        return ret

    if use_raw_data:
        return normalise_path(
            ret,
            pat_anchor_file=r"transcripts\.parquet$",
            pat_flags=re.IGNORECASE,
            return_dir=False,
            check_exist=True
        )

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_anchor_file="transcripts.parquet",
        pat_flags=re.IGNORECASE,
        return_dir=False,
        check_exist=False
    )


#######################################
#                Rules                #
#######################################

rule runProseg:
    input:
        get_input2_or_params4runProseg
    output:
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-metadata.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/transcript-metadata.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-polygons.geojson.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/expected-counts.csv.gz')
    log:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/logs/runProseg.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results',
        abs_input=lambda wildcards: os.path.abspath(
            get_input2_or_params4runProseg(
                wildcards,
                for_input=False
            )
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/logs/runProseg.log'
        ),
        other_options=get_dict_value(
            config,
            "segmentation",
            "proseg",
            "_other_options",
            replace_none=""
        )
    threads:
        get_dict_value(
            config,
            "segmentation",
            "proseg",
            "_threads"
        )
    resources:
        mem_mb=lambda wildcards, attempt: get_dict_value(
            config,
            "segmentation",
            "proseg",
            "_memory"
        ) * 1024 * attempt
    container:
        config["containers"]["proseg"]
    shell:
        "cd {params.work_dir} && "
        "proseg --nthreads {threads} "
        "{params.other_options} "
        "--xenium {params.abs_input} &> {params.abs_log}"

rule runProseg2Baysor:
    input:
        segmentation=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/transcript-metadata.csv.gz',
        polygons=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-polygons.geojson.gz'
    output:
        segmentation=protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/processed_results/baysor-transcript-metadata.csv'),
        polygons=protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/processed_results/baysor-cell-polygons.geojson')
    log:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/logs/runProseg2Baysor.log'
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 2048)
    container:
        config["containers"]["proseg"]
    shell:
        "proseg-to-baysor "
        "{input.segmentation} "
        "{input.polygons} "
        "--output-transcript-metadata {output.segmentation} "
        "--output-cell-polygons {output.polygons} &> {log}"

rule normaliseProseg:
    input:
        data_dir=lambda wildcards: get_input2_or_params4run10x(wildcards),
        segmentation=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/processed_results/baysor-transcript-metadata.csv',
        polygons=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/processed_results/baysor-cell-polygons.geojson'
    output:
        directory(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/logs/normaliseProseg.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}',
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
            f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/logs/normaliseProseg.log'
        ),
        localmem=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_memory"
        )
    retries:
        0
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
