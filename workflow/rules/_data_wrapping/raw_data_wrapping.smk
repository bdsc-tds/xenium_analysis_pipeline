#######################################
#              Functions              #
#######################################

def get_input2wrapRawDataPerSample(wildcards) -> str:
    sample_id: str = extract_layers_from_experiments(
        wildcards.wrap_sample_id,
        [0, 1, 2, 3],
        "-",
    )[0]

    return normalise_path(
        f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{sample_id}',
        pat_flags=re.IGNORECASE,
        return_dir=True,
    )


def get_input2collectWrapsAndChecksumsPerPanel(wildcards) -> dict[str, list[str]]:
    panel_id: str = extract_layers_from_experiments(
        wildcards.wrap_gene_panel_id,
        [0, 1],
        "-",
    )[0]

    sample_ids: list[str] = get_dict_value(
        config,
        "experiments",
        "_collections",
        "_gene_panels",
        panel_id,
    )

    wrap_sample_ids: list[str] = [
        extract_layers_from_experiments(
            i,
            [0, 1, 2, 3],
            sep_out="-",
        )[0]
        for i in sample_ids
    ]

    _prefix: str = f'{config["output_path"]}/wraps/raw_data/per_sample'
    wrap_suffix: str = ".tgz"
    checksum_suffix: str = ".csv"

    wrap_inputs: list[str] = [
        os.path.join(
            _prefix,
            "".join(
                [i, wrap_suffix]
            ),
        )
        for i in wrap_sample_ids
    ]
    checksum_inputs: list[str] = [
        os.path.join(
            _prefix,
            "".join(
                [i, checksum_suffix]
            ),
        )
        for i in wrap_sample_ids
    ]

    return {
        "wrap": wrap_inputs,
        "checksum": checksum_inputs,
    }


#######################################
#                Rules                #
#######################################

rule wrapRawDataPerSample:
    input:
        get_input2wrapRawDataPerSample
    output:
        f'{config["output_path"]}/wraps/raw_data/per_sample/{{wrap_sample_id}}.tgz'
    log:
        f'{config["output_path"]}/wraps/raw_data/per_sample/logs/wrap_{{wrap_sample_id}}.log'
    resources:
        runtime=300
    shell:
        "tar -czf {output} -C {input} . &> {log}"

rule computeChecksum:
    input:
        f'{config["output_path"]}/wraps/raw_data/per_sample/{{wrap_sample_id}}.tgz'
    output:
        f'{config["output_path"]}/wraps/raw_data/per_sample/{{wrap_sample_id}}.csv'
    log:
        f'{config["output_path"]}/wraps/raw_data/per_sample/logs/compute_checksum_{{wrap_sample_id}}.log'
    threads:
        1
    resources:
        runtime=60,
        mem_mb=lambda wildcards, input, attempt: min(2048 * attempt, 40960)
    shell:
        "mamba run -n utilnest python workflow/scripts/_data_wrapping/compute_checksum.py "
        "--single_file {input} "
        "-o {output} "
        "--algo sha512 "
        "--mark_algo "
        "-l {log}"

rule collectWrapsAndChecksumsPerPanel:
    input:
        unpack(get_input2collectWrapsAndChecksumsPerPanel)
    output:
        wrap=f'{config["output_path"]}/wraps/raw_data/per_gene_panel/{{wrap_gene_panel_id}}.tgz',
        checksum=f'{config["output_path"]}/wraps/raw_data/per_gene_panel/{{wrap_gene_panel_id}}.csv'
    params:
        work_dir=f'{config["output_path"]}/wraps/raw_data/per_sample'
    log:
        f'{config["output_path"]}/wraps/raw_data/per_gene_panel/logs/{{wrap_gene_panel_id}}.log'
    conda:
        "../../envs/pyarrow.yml"
    threads:
        1
    resources:
        runtime=600,
        mem_mb=lambda wildcards, attempt: min(5120 * attempt, 51200)
    script:
        "../../scripts/_data_wrapping/collect_wrapps_and_checksums_per_gene_panel.py"
