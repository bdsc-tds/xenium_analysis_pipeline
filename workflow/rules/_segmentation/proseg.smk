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

def get_input2_or_params4mapProsegRawAndNormalisedCells(wildcards, for_input: bool = True) -> list[str]:
    ret: list[str] = [
        f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/raw_results/cell-metadata.csv.gz'
    ]

    if for_input:
        ret.append(f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/normalised_results')
    else:
        ret.append(
            normalise_path(
                f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/normalised_results',
                candidate_paths=("outs",),
                pat_anchor_file="cells.parquet",
                pat_flags=re.IGNORECASE,
                return_dir=False,
                check_exist=False
            )
        )

    return ret


#######################################
#                Rules                #
#######################################

rule runProseg:
    input:
        get_input2_or_params4runProseg
    output:
        directory(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/proseg-output.zarr'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/counts.mtx.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/expected-counts.mtx.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-metadata.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/transcript-metadata.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/gene-metadata.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-polygons.geojson.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-polygons-layers.geojson.gz')
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
            "run",
            "_other_options",
            replace_none="",
        )
    threads:
        get_dict_value(
            config,
            "segmentation",
            "proseg",
            "run",
            "_threads",
            replace_none=1,
        )
    resources:
        mem_mb=lambda wildcards, attempt: get_dict_value(
            config,
            "segmentation",
            "proseg",
            "run",
            "_memory",
            replace_none=20,
        ) * 1024 * attempt
    container:
        config["containers"]["proseg"]
    shell:
        "cd {params.work_dir} && "
        "proseg --nthreads {threads} "
        "{params.other_options} "
        "--output-counts counts.mtx.gz "
        "--output-expected-counts expected-counts.mtx.gz "
        "--output-cell-metadata cell-metadata.csv.gz "
        "--output-transcript-metadata transcript-metadata.csv.gz "
        "--output-gene-metadata gene-metadata.csv.gz "
        "--output-cell-polygons cell-polygons.geojson.gz "
        "--output-cell-polygon-layers cell-polygons-layers.geojson.gz "
        "--xenium {params.abs_input} &> {params.abs_log}"

rule runProseg2Baysor:
    input:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/proseg-output.zarr'
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
        "{input} "
        "--output-transcript-metadata {output.segmentation} "
        "--output-cell-polygons {output.polygons} &> {log}"

rule convertProsegCountsFormat:
    input:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/counts.mtx.gz',
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/expected-counts.mtx.gz',
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/gene-metadata.csv.gz'
    output:
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/counts.csv.gz'),
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/expected-counts.csv.gz')
    log:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/logs/runProseg.log'
    threads:
        1
    conda:
        "../../envs/coexpression.yml"
    shell:
        "python3 workflow/scripts/_segmentation/proseg_mtx_to_csv.py "
        "--counts {input[0]} "
        "--expected_counts {input[1]} "
        "--gene {input[2]} "
        "--out_counts {output[0]} "
        "--out_expected_counts {output[1]} "
        "-l {log}"

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

rule mapProsegRawAndNormalisedCells:
    input:
        lambda wildcards: get_input2_or_params4mapProsegRawAndNormalisedCells(
            wildcards,
            for_input=True,
        )
    output:
        protected(f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/mapped_cell_ids/mapped_cell_ids.parquet')
    params:
        lambda wildcards: get_input2_or_params4mapProsegRawAndNormalisedCells(
            wildcards,
            for_input=False,
        )
    log:
        f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/logs/mapProsegRawAndNormalisedCells.log'
    resources:
        mem_mb=lambda wildcards, attempt: max(
            sum(
                get_size(i) for i in get_input2_or_params4mapProsegRawAndNormalisedCells(
                    wildcards,
                    for_input=True,
                )
            ) * 10**-6 * attempt * 10,
            2048
        )
    conda:
        "../../envs/pyarrow.yml"
    shell:
        "python3 workflow/scripts/_segmentation/map_proseg_raw_and_normalised_cells.py "
        "--raw {params[0][0]} "
        "--normalised {params[0][1]} "
        "-l {log} "
        "--out {output}"
