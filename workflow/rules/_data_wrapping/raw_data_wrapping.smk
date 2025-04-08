#######################################
#              Functions              #
#######################################

def get_input2wrapRawData(wildcards) -> str:
    sample_id: str = extract_layers_from_experiments(
        wildcards.wrap_sample_id,
        [0, 1, 2, 3],
        "+"
    )[0]

    return normalise_path(
        f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{sample_id}',
        pat_flags=re.IGNORECASE,
        return_dir=True,
    )


#######################################
#                Rules                #
#######################################

rule wrapRawData:
    input:
        get_input2wrapRawData
    output:
        f'{config["output_path"]}/wraps/raw_data/{{wrap_sample_id}}.tgz'
    log:
        f'{config["output_path"]}/wraps/raw_data/logs/{{wrap_sample_id}}.log'
    resources:
        runtime=300
    shell:
        "tar -czf {output} -C {input} . &> {log}"
