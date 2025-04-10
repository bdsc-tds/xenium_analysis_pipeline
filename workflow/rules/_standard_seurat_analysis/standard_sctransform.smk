#######################################
#                Rules                #
#######################################

"""
Comment:

- future_globals_maxSize: By default 10GB is used to deal with the [issue](https://github.com/satijalab/seurat/issues/1845) of "Global size exceeds maximum allowed size" when running Seurat (see the solution [here](https://satijalab.org/seurat/archive/v3.0/future_vignette.html)).
"""
rule runStandardScTransform:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/normalised_counts/scale_data.parquet')
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        assay=sec.SEURAT_ALT_ASSAY,
        data_layer=sec.SEURAT_DATA_LAYER,
        scale_data_layer=sec.SEURAT_SCALE_DATA_LAYER,
        normalisation_id="sctransform"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(input.size_mb * attempt**3 * 40, 1024000),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/logs/runStandardScTransform.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_sctransform.R"

use rule runStandardScTransform as runPostCountCorrectionStandardScTransform with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/sctransform/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/sctransform/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/sctransform/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/sctransform/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/{{annotation_id}}/{{count_correction_id}}/sctransform/logs/runPostCountCorrectionStandardScTransform.log'
    wildcard_constraints:
        count_correction_id=COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT

use rule runStandardScTransform as runPostCountCorrectionByOvrlpyStandardScTransform with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/sctransform/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/sctransform/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/sctransform/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/sctransform/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/sctransform/logs/runPostCountCorrectionByOvrlpyStandardScTransform.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule runStandardScTransform as runPostCountCorrectionByResolviUnsupervisedStandardScTransform with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/logs/runPostCountCorrectionByResolviUnsupervisedStandardScTransform.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised"

use rule runStandardScTransform as runPostCountCorrectionByResolviSupervisedStandardScTransform with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/normalised_seurat.rds',
        cells=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/cells.parquet'),
        data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/data.parquet'),
        scale_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/normalised_counts/scale_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/sctransform/logs/runPostCountCorrectionByResolviSupervisedStandardScTransform.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
