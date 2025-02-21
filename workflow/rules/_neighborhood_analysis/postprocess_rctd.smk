#######################################
#              Functions              #
#######################################



#######################################
#                Rules                #
#######################################

rule runPostprocessRCTD:
    input:
        rctd_result=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'
    output:
        post_processed_rctd=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/post_processed_output.rds'),
        post_processed_rctd_df=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/post_processed_results_df.parquet')
    params:
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedRCTDpostprocessing.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_neighborhood_analysis/postprocess_rctd.R"
