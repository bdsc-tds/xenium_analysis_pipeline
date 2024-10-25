#######################################
#                Rules                #
#######################################

rule generateStandardSeuratReport:
    input:
        raw=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/raw_seurat.rds',
        preprocessed=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/preprocessed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/reports/standard_seurat_analysis.html')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        sample_id=lambda wildcards: wildcards.sample_id,
        rmd_file="workflow/scripts/_standard_seurat_analysis/standard_seurat_report.Rmd"
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/generateStandardSeuratReport.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/generate_standard_seurat_report.R"
