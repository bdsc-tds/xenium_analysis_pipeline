#######################################
#                Rules                #
#######################################

rule computeCoexpression:
    input:
        get_input2_or_params4loadSegmentation2Seurat
    output:
        co_exp=protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/coexpression/{{coexpression_id}}/coexpression.parquet'),
        pos_rate=protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/coexpression/{{coexpression_id}}/positivity_rate.parquet')
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/coexpression/{{coexpression_id}}/logs/computeCoexpression.log'
    params:
        data_dir=lambda wildcards: get_input2_or_params4loadSegmentation2Seurat(
            wildcards,
            for_input=False
        ),
        method=lambda wildcards: extract_layers_from_experiments(
            wildcards.coexpression_id,
            layers=[0],
            sep_in="_",
        )[0],
        target_counts=lambda wildcards: extract_layers_from_experiments(
            wildcards.coexpression_id,
            layers=[1],
            sep_in="_",
        )[0]
    container:
        config["containers"]["r"]
    conda:
        "../envs/coexpression.yml"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt**4 * 512
    shell:
        "python3 workflow/scripts/_coexpression/compute_coexpression.py "
        "-i {params.data_dir} "
        "-l {log} "
        "--outcoexp {output.co_exp} "
        "--outposrate {output.pos_rate} "
        "-m {params.method} "
        "-c {params.target_counts}"
