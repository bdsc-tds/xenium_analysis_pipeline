#######################################
#                Rules                #
#######################################

rule runStandardQC:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds'
    output:
        obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds',
        meta_data=protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/meta_data.parquet')
    params:
        min_counts=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_counts",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_counts",
                replace_none=10
            ),
            inexist_key_ok=True
        ),
        min_features=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_features",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_features",
                replace_none=5
            ),
            inexist_key_ok=True
        ),
        max_counts=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "max_counts",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "max_counts",
                replace_none="Inf"
            ),
            inexist_key_ok=True
        ),
        max_features=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "max_features",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "max_features",
                replace_none="Inf"
            ),
            inexist_key_ok=True
        ),
        min_cells=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0, 1])[0],
            "min_cells",
            replace_none=get_dict_value(
                config,
                "standard_seurat_analysis",
                "qc",
                "min_cells",
                replace_none=5
            ),
            inexist_key_ok=True
        ),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        default_layer=sec.SEURAT_DEFAULT_LAYER
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/runStandardQC.log'
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 20480)
    script:
        "../../scripts/_standard_seurat_analysis/standard_qc.R"

use rule runStandardQC as runPostCountCorrectionStandardQC with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/raw_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/standard_qc/qced_seurat.rds',
        meta_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/standard_qc/meta_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/logs/runPostCountCorrectionStandardQC.log'
    wildcard_constraints:
        count_correction_id=COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT

use rule runStandardQC as runPostCountCorrectionByOvrlpyStandardQC with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/raw_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/standard_qc/qced_seurat.rds',
        meta_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/standard_qc/meta_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/logs/runPostCountCorrectionByOvrlpyStandardQC.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule runStandardQC as runPostCountCorrectionByResolviUnsupervisedStandardQC with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/raw_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds',
        meta_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/meta_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/logs/runPostCountCorrectionByResolviUnsupervisedStandardQC.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised"

use rule runStandardQC as runPostCountCorrectionByResolviSupervisedStandardQC with:
    input:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/raw_seurat.rds'
    output:
        obj=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/qced_seurat.rds',
        meta_data=protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/standard_qc/meta_data.parquet')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/logs/runPostCountCorrectionByResolviSupervisedStandardQC.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
