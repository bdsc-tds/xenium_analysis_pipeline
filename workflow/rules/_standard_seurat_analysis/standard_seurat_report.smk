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
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 10240)
    log:
        f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/logs/generateStandardSeuratReport.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/generate_standard_seurat_report.R"

use rule generateStandardSeuratReport as generatePostCountCorrectionStandardSeuratReport with:
    input:
        raw=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/raw_seurat.rds',
        preprocessed=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/standard_seurat_analysis.html')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        sample_id=lambda wildcards: wildcards.sample_id,
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        rmd_file="workflow/scripts/_standard_seurat_analysis/standard_seurat_report.Rmd",
        intermediates_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/_intermediates_seurat_report',
        knit_root_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/_knit_root_seurat_report'
    log:
        f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/logs/generatePostCountCorrectionStandardSeuratReport.log'
    wildcard_constraints:
        count_correction_id=COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT

use rule generateStandardSeuratReport as generatePostCountCorrectionByOvrlpyStandardSeuratReport with:
    input:
        raw=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/raw_seurat.rds',
        preprocessed=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/standard_seurat_analysis.html')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        sample_id=lambda wildcards: wildcards.sample_id,
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        rmd_file="workflow/scripts/_standard_seurat_analysis/standard_seurat_report.Rmd",
        intermediates_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/_intermediates_seurat_report',
        knit_root_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/_knit_root_seurat_report'
    log:
        f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/logs/generatePostCountCorrectionByOvrlpyStandardSeuratReport.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule generateStandardSeuratReport as generatePostCountCorrectionByResolviUnsupervisedStandardSeuratReport with:
    input:
        raw=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/raw_seurat.rds',
        preprocessed=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'
    output:
        protected(f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/standard_seurat_analysis.html')
    params:
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        sample_id=lambda wildcards: wildcards.sample_id,
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        rmd_file="workflow/scripts/_standard_seurat_analysis/standard_seurat_report.Rmd",
        intermediates_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/_intermediates_seurat_report',
        knit_root_dir=f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/_knit_root_seurat_report'
    log:
        f'{config["output_path"]}/reports/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/{{normalisation_id}}/logs/generatePostCountCorrectionByResolviUnsupervisedStandardSeuratReport.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised"
