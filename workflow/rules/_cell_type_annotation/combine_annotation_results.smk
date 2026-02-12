#######################################
#                Rules                #
#######################################

rule combineAnnotationLabels:
    input:
        # Directory to search for tool/level/mode label files
        root_dir=lambda wildcards: os.path.join(
            config["output_path"],
            "cell_type_annotation",
            wildcards.segmentation_id,
            wildcards.sample_id,
            wildcards.normalisation_id,
            "reference_based",
            wildcards.reference_name
        )
    output:
        combined=protected(
            f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/reference_based/{{reference_name}}/combined_labels.parquet'
        )
    params:
        script="scripts/_cell_type_annotation/combine_annotation_results.R"
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/reference_based/{{reference_name}}/logs/combineAnnotationLabels.log'
    container:
        config["containers"]["r"]
    threads: 1
    shell:
        """
        Rscript {params.script} \
            --input_dir {input.root_dir} \
            --output_file {output.combined} \
            &> {log}
        """
