#######################################
#                Rules                #
#######################################

rule generateSegmentationQCReport:
    input:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/segmentation_qc.csv'
    output:
        protected(f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/segmentation_qc.html')
    params:
        rmd_file = "workflow/scripts/_segmentation_qc/segmentation_qc_report.Rmd"
    log:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel}}/logs/generateSegmentationQCReport.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_segmentation_qc/generate_segmentation_qc_report.R"
