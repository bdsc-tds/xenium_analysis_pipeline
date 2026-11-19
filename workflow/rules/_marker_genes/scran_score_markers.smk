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
        min_cells_per_type=get_marker_genes_config(
            "min_cells_per_type",
            default=25,
        ),
        rctd_singlets_only=get_marker_genes_config(
            "rctd_singlets_only",
            default=True,
        )
    wildcard_constraints:
        marker_id=r"scran_score_markers"
    threads:
        get_marker_genes_config(
            "scran_score_markers",
            "_threads",
            default=4,
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
