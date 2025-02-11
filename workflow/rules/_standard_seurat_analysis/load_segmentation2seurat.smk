#######################################
#              Functions              #
#######################################

def get_input2_or_params4loadSegmentation2Seurat(wildcards, for_input: bool = True) -> str:
    seg_method: str = "proseg" if wildcards.segmentation_id.startswith("proseg") else wildcards.segmentation_id

    ret = f'{config["output_path"]}/segmentation/{seg_method}/{wildcards.sample_id}/normalised_results'

    if for_input:
        return ret

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False
    )

def use_mode_counts4loadProseg2Seurat(wildcards) -> bool:
    matched = re.match(r"^proseg_(\w+)", wildcards.segmentation_id, flags=re.IGNORECASE)
    assert matched is not None

    if matched.group(1) == "mode":
        return True
    elif matched.group(1) == "expected":
        return False

    raise RunError(f"Error! Invalid segmentation method: {matched.group(1)}")

def get_mapping_param4loadProseg2Seurat(wildcards) -> bool:
    use_mode_counts: bool = use_mode_counts4loadProseg2Seurat(wildcards)

    return get_dict_value(
        config,
        "segmentation",
        "proseg",
        "mode" if use_mode_counts else "expected",
        "use_mapping",
    )


#######################################
#                Rules                #
#######################################

rule loadSegmentation2Seurat:
    input:
        get_input2_or_params4loadSegmentation2Seurat
    output:
        protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds')
    params:
        data_dir=lambda wildcards: get_input2_or_params4loadSegmentation2Seurat(
            wildcards,
            for_input=False
        ),
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        sample_id=lambda wildcards: wildcards.sample_id,
        segmentation_id=lambda wildcards: wildcards.segmentation_id
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadSegmentation2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        segmentation_id=r"(?!proseg).*"
    script:
        "../../scripts/_standard_seurat_analysis/load_segmentation2seurat.R"

rule loadProseg2Seurat:
    input:
        dir=get_input2_or_params4loadSegmentation2Seurat,
        mapping=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/mapped_cell_ids/mapped_cell_ids.parquet',
        cell_metadata=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/cell-metadata.csv.gz',
        expected_counts=f'{config["output_path"]}/segmentation/proseg/{{sample_id}}/raw_results/expected-counts.csv.gz'
    output:
        protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds')
    params:
        data_dir=lambda wildcards: get_input2_or_params4loadSegmentation2Seurat(
            wildcards,
            for_input=False
        ),
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        sample_id=lambda wildcards: wildcards.sample_id,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        xr_cell_id_col_name="xr_cell_id",
        proseg_cell_id_col_name="proseg_cell_id",
        use_mode_counts=use_mode_counts4loadProseg2Seurat,
        use_mapping=get_mapping_param4loadProseg2Seurat
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadProseg2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        segmentation_id=r"proseg_\w+"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(
            sum(
                get_size(i) for i in input if os.path.isfile(i)
            ) * 10**-6 * attempt * 4,
            2048
        )
    script:
        "../../scripts/_standard_seurat_analysis/load_proseg2seurat.R"
