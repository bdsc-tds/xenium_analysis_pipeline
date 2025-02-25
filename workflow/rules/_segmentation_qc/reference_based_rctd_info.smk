#######################################
#              Functions              #
#######################################

def generate_dummy_table(wildcards) -> bool:
    files: list[str] = [
        f'{config["output_path"]}/std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}/standard_qc/meta_data.parquet',
        f'{config["output_path"]}/cell_type_annotation/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}/{wildcards.annotation_id}/output/results_df.parquet',
    ]

    if any(
        not os.path.exists(file)
        for file in files
    ) and get_dict_value(
        config,
        "segmentation",
        "_qc",
        "_create_dummy_table",
        replace_none=False,
    ):
        return True
    
    return False


def get_input2gatherReferenceBasedRCTDInfo(wildcards) -> str | list[str]:
    if generate_dummy_table(wildcards):
        return f'{config["output_path"]}/reprocessed/{wildcards.sample_id}/versions.json'

    return [
        f'{config["output_path"]}/std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}/standard_qc/meta_data.parquet',
        f'{config["output_path"]}/cell_type_annotation/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}/{wildcards.annotation_id}/output/results_df.parquet',
    ]


def get_shell_cmd4gatherReferenceBasedRCTDInfo(
    wildcards,
    input,
    output,
    log,
) -> str:
    segmentation_method = extract_layers_from_experiments(
        wildcards.segmentation_id,
        0,
        sep_in="_",
    )[0]

    condition = extract_layers_from_experiments(
        wildcards.sample_id,
        0,
    )[0]
    gene_panel = extract_layers_from_experiments(
        wildcards.sample_id,
        1,
    )[0]
    donor = extract_layers_from_experiments(
        wildcards.sample_id,
        2,
    )[0]
    sample = extract_layers_from_experiments(
        wildcards.sample_id,
        3,
    )[0]

    reference_name = extract_layers_from_experiments(
        wildcards.annotation_id,
        1,
    )[0]
    reference_level = extract_layers_from_experiments(
        wildcards.annotation_id,
        3,
    )[0]

    if generate_dummy_table(wildcards):
        return " ".join(
            [
                "python3",
                "workflow/scripts/_segmentation_qc/create_dummy_table.py",
                "--segmentation_id",
                wildcards.segmentation_id,
                "--segmentation_method",
                segmentation_method,
                "--sample_id",
                wildcards.sample_id,
                "--condition",
                condition,
                "--gene_panel",
                gene_panel,
                "--donor",
                donor,
                "--sample",
                sample,
                "--normalisation_id",
                wildcards.normalisation_id,
                "--annotation_id",
                wildcards.annotation_id,
                "--reference_name",
                reference_name,
                "--reference_level",
                reference_level,
                "--na_values",
                "--out",
                output[0],
                "-l",
                log,
            ]
        )

    assert isinstance(input, list) and len(input) == 2
    return " ".join(
        [
            "python3",
            "workflow/scripts/_segmentation_qc/gather_reference_based_rctd_info.py",
            "--in_meta",
            input[0],
            "--in_rctd",
            input[1],
            "--out",
            output[0],
            "-l",
            log,
            "--normalisation_id",
            wildcards.normalisation_id,
            "--annotation_id",
            wildcards.annotation_id,
            "--reference_name",
            reference_name,
            "--reference_level",
            reference_level,
        ]
    )


#######################################
#                Rules                #
#######################################

rule gatherReferenceBasedRCTDInfo:
    input:
        get_input2gatherReferenceBasedRCTDInfo
    output:
        protected(f'{config["output_path"]}/segmentation_qc/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_info.parquet')
    params:
        cmd=lambda wildcards, input, output: get_shell_cmd4gatherReferenceBasedRCTDInfo(
            wildcards,
            input,
            output,
            f'{config["output_path"]}/segmentation_qc/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}/{wildcards.annotation_id}/logs/gatherReferenceBasedRCTDInfo.log',
        )
    log:
        f'{config["output_path"]}/segmentation_qc/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/gatherReferenceBasedRCTDInfo.log'
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd_.+"
    conda:
        "../../envs/pyarrow.yml"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(
            input.size_mb * attempt * 100,
            256,
        )
    shell:
        "{params.cmd}"
