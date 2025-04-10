#######################################
#                Rules                #
#######################################

rule runStandardLogNorm:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/normalised_counts/scale_data.parquet')
    params:
        assay=sec.SEURAT_DEFAULT_ASSAY,
        data_layer=sec.SEURAT_DATA_LAYER,
        scale_data_layer=sec.SEURAT_SCALE_DATA_LAYER,
        normalisation_id="lognorm"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt**3 * 40, 10240)
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/logs/runStandardLogNorm.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_lognorm.R"

use rule runStandardLogNorm as runPostCountCorrectionStandardLogNorm with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/lognorm/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/lognorm/{{annotation_id}}/{{count_correction_id}}/lognorm/logs/runPostCountCorrectionStandardLogNorm.log'
    wildcard_constraints:
        count_correction_id=COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT

use rule runStandardLogNorm as runPostCountCorrectionByOvrlpyStandardLogNorm with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/lognorm/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/lognorm/logs/runPostCountCorrectionByOvrlpyStandardLogNorm.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule runStandardLogNorm as runPostCountCorrectionByResolviUnsupervisedStandardLogNorm with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/logs/runPostCountCorrectionByResolviUnsupervisedStandardLogNorm.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised"

use rule runStandardLogNorm as runPostCountCorrectionByResolviSupervisedStandardLogNorm with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/lognorm/logs/runPostCountCorrectionByResolviUnsupervisedStandardLogNorm.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
