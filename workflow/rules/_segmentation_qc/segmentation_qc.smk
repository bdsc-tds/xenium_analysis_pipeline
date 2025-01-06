#######################################
#                Rules                #
#######################################

rule segmentationQC:
    input:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/metadata.csv.gz'
    output:
        protected(f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/segmentation_qc.csv')
    params:
    log:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/logs/segmentationQC.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_segmentation_qc/segmentation_qc.R"
