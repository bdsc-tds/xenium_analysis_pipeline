#######################################
#                Rules                #
#######################################

rule generateStandardSeuratReport:
    input:
        raw=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds',
        preprocessed=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/standard_seurat_analysis.html')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        sample_id=lambda wildcards: wildcards.sample_id,
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        rmd_file="workflow/scripts/_standard_seurat_analysis/standard_seurat_report.Rmd",
        intermediates_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/_intermediates_seurat_report',
        knit_root_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/_knit_root_seurat_report'
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 30, 10240)
    log:
        f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/logs/generateStandardSeuratReport.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/generate_standard_seurat_report.R"
