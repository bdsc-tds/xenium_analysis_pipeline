#######################################
#                Rules                #
#######################################

rule runScranScoreMarkers:
    input:
        xe=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        labels=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        rctd_results=lambda wildcards: get_rctd_results4marker_genes(
            wildcards,
            f'{config["output_path"]}/cell_type_annotation',
        )
    output:
        markers=f'{config["output_path"]}/marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/markers.parquet',
        pseudobulk=f'{config["output_path"]}/marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/pseudobulk.parquet'
    params:
        xe_assay=lambda wildcards: get_assay_name4annotation(wildcards, False),
        counts_assay=sec.SEURAT_DEFAULT_ASSAY,
        control_gene_pat=sec.XENIUM_CONTROL_GENE_PAT,
        min_cells_per_type=get_dict_value(
            config,
            "marker_genes",
            "min_cells_per_type",
            replace_none=25,
            inexist_key_ok=True,
        ),
        rctd_singlets_only=get_dict_value(
            config,
            "marker_genes",
            "rctd_singlets_only",
            replace_none=True,
            inexist_key_ok=True,
        )
    wildcard_constraints:
        marker_id=r"scran_score_markers"
    threads:
        get_dict_value(
            config,
            "marker_genes",
            "scran_score_markers",
            "_threads",
            replace_none=4,
            inexist_key_ok=True,
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            max(input.size_mb * attempt * 20, 10240),
            1024000,
        )
    log:
        f'{config["output_path"]}/marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/logs/runScranScoreMarkers.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_marker_genes/scran_score_markers.R"

use rule runScranScoreMarkers as runPostCountCorrectionBySplitScranScoreMarkers with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        labels=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        rctd_results=lambda wildcards: get_rctd_results4marker_genes(
            wildcards,
            f'{config["output_path"]}/post_count_correction_cell_type_annotation',
        )
    output:
        markers=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/markers.parquet',
        pseudobulk=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/pseudobulk.parquet'
    log:
        f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/logs/runPostCountCorrectionBySplitScranScoreMarkers.log'
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r".+/rctd_.+",
        marker_id=r"scran_score_markers"

use rule runScranScoreMarkers as runPostCountCorrectionByOvrlpyScranScoreMarkers with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        labels=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        rctd_results=lambda wildcards: get_rctd_results4marker_genes(
            wildcards,
            f'{config["output_path"]}/post_count_correction_cell_type_annotation',
        )
    output:
        markers=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/markers.parquet',
        pseudobulk=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/pseudobulk.parquet'
    log:
        f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/logs/runPostCountCorrectionByOvrlpyScranScoreMarkers.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        marker_id=r"scran_score_markers"

use rule runScranScoreMarkers as runPostCountCorrectionByResolviUnsupervisedScranScoreMarkers with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        labels=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        rctd_results=lambda wildcards: get_rctd_results4marker_genes(
            wildcards,
            f'{config["output_path"]}/post_count_correction_cell_type_annotation',
        )
    output:
        markers=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/markers.parquet',
        pseudobulk=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/pseudobulk.parquet'
    log:
        f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/logs/runPostCountCorrectionByResolviUnsupervisedScranScoreMarkers.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        marker_id=r"scran_score_markers"

use rule runScranScoreMarkers as runPostCountCorrectionByResolviSupervisedScranScoreMarkers with:
    input:
        xe=f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds',
        labels=f'{config["output_path"]}/post_count_correction_cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        rctd_results=lambda wildcards: get_rctd_results4marker_genes(
            wildcards,
            f'{config["output_path"]}/post_count_correction_cell_type_annotation',
        )
    output:
        markers=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/markers.parquet',
        pseudobulk=f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/pseudobulk.parquet'
    log:
        f'{config["output_path"]}/post_count_correction_marker_genes/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/{{normalisation_id}}/{{annotation_id}}/{{marker_id}}/logs/runPostCountCorrectionByResolviSupervisedScranScoreMarkers.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        marker_id=r"scran_score_markers"
