#######################################
#              Functions              #
#######################################

def get_input2check10xVersions(wildcards) -> str:
    return normalise_path(
        f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{wildcards.sample_id}',
        pat_flags=re.IGNORECASE,
        return_dir=False
    )

def get_gene_panel4reprocessRawData(wildcards) -> str:
    ret: str | None = get_gene_panel_file(wildcards.sample_id, config)

    if ret is None:
        raise RuntimeError(f"Error! Sample {wildcards.sample_id} does not have related gene panel file provided.")

    return ret

def get_input2_or_params4changeParquetCompressionType(wildcards, for_input: bool = True) -> str:
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

checkpoint check10xVersions:
    input:
        get_input2check10xVersions
    output:
        protected(f'{config["output_path"]}/reprocessed/{{sample_id}}/versions.json')
    log:
        f'{config["output_path"]}/reprocessed/{{sample_id}}/logs/check10xVersions.log'
    resources:
        mem_mb=1024
    conda:
        "../envs/bs4.yml"
    container:
        config["containers"]["10x"]
    shell:
        "python3 workflow/scripts/collect_10x_versions.py "
        "-i {input} "
        "-l {log} "
        "-o {output}"

rule reprocessRawData:
    input:
        data_dir=f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{{sample_id}}'
    output:
        directory(f'{config["output_path"]}/reprocessed/{{sample_id}}/results')
    log:
        f'{config["output_path"]}/reprocessed/{{sample_id}}/logs/reprocessRawData.log'
    params:
        work_dir=f'{config["output_path"]}/reprocessed/{{sample_id}}',
        abs_input=lambda wildcards, input: os.path.abspath(
            normalise_path(
                input["data_dir"],
                pat_flags=re.IGNORECASE
            )
        ),
        abs_gene_panel=lambda wildcards: get_gene_panel4reprocessRawData(wildcards),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/reprocessed/{wildcards.sample_id}/logs/reprocessRawData.log'
        ),
        localmem=lambda wildcards: get_dict_value(
            config,
            "reprocess",
            "_memory"
        )
    threads:
        lambda wildcards: get_dict_value(
            config,
            "reprocess",
            "_threads"
        )
    resources:
        mem_mb=lambda wildcards: get_dict_value(
            config,
            "reprocess",
            "_memory"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger relabel --id=results "
        "--xenium-bundle={params.abs_input} "
        "--panel={params.abs_gene_panel} "
        "--localcores={threads} "
        "--localmem={params.localmem} &> {params.abs_log}"

rule changeParquetCompressionType:
    input:
        get_input2_or_params4changeParquetCompressionType
    output:
        protected(f'{config["output_path"]}/reprocessed/{{sample_id}}/transcripts_snappy.parquet')
    log:
        f'{config["output_path"]}/reprocessed/{{sample_id}}/logs/changeParquetCompressionType.log'
    params:
        input_file=lambda wildcards: get_input2_or_params4changeParquetCompressionType(wildcards, for_input=False)
    resources:
        mem_mb=lambda wildcards, input: max(input.size_mb * 10, 20480)
    conda:
        "../envs/pyarrow.yml"
    shell:
        "python3 workflow/scripts/change_parquet_compresson_type.py "
        "-i {params.input_file} "
        "-t snappy "
        "-l {log} "
        "-o {output}"
