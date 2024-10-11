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
