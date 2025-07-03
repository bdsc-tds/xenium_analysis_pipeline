#######################################
#              Functions              #
#######################################

def get_input2_or_params4runBats(wildcards, for_input: bool = True) -> str:
    seg_method: atr = get_dict_value(
        config,
        "segmentation",
        "bats",
        "load_from",
    )
    
    ret = f'{config["output_path"]}/segmentation/{seg_method}/{wildcards.sample_id}/normalised_results'

    if for_input:
        return ret

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False
    )


#######################################
#                Rules                #
#######################################

rule runBats:
    input:
        get_input2_or_params4runBats
    output:
        protected(f'{config["output_path"]}/segmentation/bats/{{sample_id}}/raw_results/expected_counts.parquet')
    log:
        f'{config["output_path"]}/segmentation/bats/{{sample_id}}/logs/runBats.log'
    params:
        input=lambda wildcards: get_input2_or_params4runBats(
            wildcards,
            for_input=False,
        ),
        other_options=get_dict_value(
            config,
            "segmentation",
            "bats",
            "_other_options",
            replace_none="",
        )
    threads:
        get_dict_value(
            config,
            "segmentation",
            "bats",
            "_threads",
            replace_none=1,
        )
    resources:
        mem_mb=lambda wildcards, attempt: get_dict_value(
            config,
            "segmentation",
            "bats",
            "_memory",
            replace_none=2,
        ) * 1024 * attempt
    # container:
    #     config["containers"]["python_cuda"]
    shell:
        "echo 'Simulating running Bats...' > {log}"
