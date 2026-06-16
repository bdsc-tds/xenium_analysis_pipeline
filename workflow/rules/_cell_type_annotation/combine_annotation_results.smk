#######################################
#                Rules                #
#######################################

rule combineAnnotationLabels:
    wildcard_constraints:
        reference_name="[^/]+"
        
    input:
      root_dir=lambda wc: os.path.join(
            config["output_path"],
            "cell_type_annotation",
            wc.segmentation_id,
            wc.sample_id,
            wc.normalisation_id,
            "reference_based",
            wc.reference_name
        )

    output:
        combined=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/reference_based/{{reference_name}}/combined_labels.parquet'
    
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/reference_based/{{reference_name}}/logs/combineAnnotationLabels.log'
    container:
        config["containers"]["r"]
    threads: 1
    script:
        "../../scripts/_cell_type_annotation/combine_annotation_results.R"

