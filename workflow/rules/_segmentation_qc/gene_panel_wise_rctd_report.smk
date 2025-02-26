#######################################
#                Rules                #
#######################################

rule generateGenePanelWiseRCTDReport:
    input:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel_id}}/rctd_info.parquet'
    output:
        protected(f'{config["output_path"]}/reports/{{gene_panel_id}}/segmentation_qc_rctd.html')
    params:
        gene_panel_id=lambda wildcards: wildcards.gene_panel_id,
        rmd_file="workflow/scripts/_segmentation_qc/segmentation_qc_rctd_report.Rmd",
        reference_name=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "_qc",
            "reference_name",
        ),
        reference_level=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "_qc",
            "reference_level",
        ),
        intermediates_dir=f'{config["output_path"]}/reports/{{gene_panel_id}}/_intermediates_rctd_report',
        knit_root_dir=f'{config["output_path"]}/reports/{{gene_panel_id}}/_knit_root_rctd_report'
    log:
        f'{config["output_path"]}/reports/{{gene_panel_id}}/logs/generateGenePanelWiseRCTDReport.log'
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 1024, 20480)
    script:
        "../../scripts/_segmentation_qc/generate_rctd_report.R"
