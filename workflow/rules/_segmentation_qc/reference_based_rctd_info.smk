#######################################
#                Rules                #
#######################################

rule gatherReferenceBasedRCTDInfo:
    input:
        sample_info=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/meta_data.parquet',
        anno_results=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet'
    output:
        protected(f'{config["output_path"]}/segmentation_qc/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_info.parquet')
    params:
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        annotation_id=lambda wildcards: wildcards.annotation_id,
        ref_name=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            1,
        )[0],
        ref_level=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            3,
        )[0]
    log:
        f'{config["output_path"]}/segmentation_qc/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/gatherReferenceBasedRCTDInfo.log'
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd_.+"
    conda:
        "../../envs/pyarrow.yml"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(
            input.size_mb * attempt * 100,
            20480,
        )
    shell:
        "python3 workflow/scripts/_segmentation_qc/gather_reference_based_rctd_info.py "
        "--in_meta {input.sample_info} "
        "--in_rctd {input.anno_results} "
        "--out {output} "
        "-l {log} "
        "--normalisation_id {params.normalisation_id} "
        "--annotation_id {params.annotation_id} "
        "--ref_name {params.ref_name} "
        "--ref_level {params.ref_level}"
