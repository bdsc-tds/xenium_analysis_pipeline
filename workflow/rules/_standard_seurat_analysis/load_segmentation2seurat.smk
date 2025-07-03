#######################################
#              Functions              #
#######################################

def get_input2_or_params4loadSegmentation2Seurat(wildcards, for_input: bool = True) -> str:
    if wildcards.segmentation_id.startswith("proseg"):
        seg_method = "proseg"
    elif wildcards.segmentation_id == "bats":
        seg_method = get_dict_value(
            config,
            "segmentation",
            "bats",
            "load_from",
        )
    else:
        seg_method = wildcards.segmentation_id

    ret = f'{config["output_path"]}/segmentation/{seg_method}/{wildcards.sample_id}/normalised_results'

    if for_input:
        return ret

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False,
    )

def get_mapped_cell_ids(
    sample_id: str,
    segmentation_id: str,
    normalised_path: str,
) -> tuple[bool, str]:
    is_xr_v4: bool = validate_xr_version(
        sample_id,
        True,
        "0",
        lambda x: x >= 4,
    )

    if is_xr_v4:
        ret:str = normalise_path(
            normalised_path,
            candidate_paths=(".",),
            pat_anchor_file="cell_id_map.csv.gz",
            pat_flags=re.IGNORECASE,
            return_dir=False,
            check_exist=False,
        )
    else:
        ret = f'{config["output_path"]}/segmentation/{segmentation_id}/{sample_id}/mapped_cell_ids/mapped_cell_ids.parquet'
    
    return is_xr_v4, ret

def get_input2_or_params4loadProseg2Seurat(wildcards, for_input: bool = True) -> dict[str, str]:
    ret: dict[str, str] = {
        "cell_metadata": f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/raw_results/cell-metadata.csv.gz',
        "expected_counts": f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/raw_results/expected-counts.csv.gz',
        "dir": get_input2_or_params4loadSegmentation2Seurat(
            wildcards,
            for_input=for_input,
        ),
    }

    is_xr_v4, path2mapping = get_mapped_cell_ids(
        wildcards.sample_id,
        "proseg",
        ret["dir"],
    )

    if not (is_xr_v4 and for_input):
        ret["mapping"] = path2mapping

    return ret

def use_mode_counts4loadProseg2Seurat(wildcards) -> bool:
    matched = re.match(
        r"^proseg_(\w+)",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    )
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

def get_imported_cell_id4mapping(wildcards) -> str:
    seg_method: str = "proseg" if wildcards.segmentation_id.startswith("proseg") else wildcards.segmentation_id

    is_xr_v4: bool = validate_xr_version(
        wildcards.sample_id,
        True,
        "0",
        lambda x: x >= 4,
    )

    return "imported_cell_id" if is_xr_v4 else "imported_cell_id"

def get_xenium_cell_id4mapping(wildcards) -> str:
    seg_method: str = "proseg" if wildcards.segmentation_id.startswith("proseg") else wildcards.segmentation_id

    is_xr_v4: bool = validate_xr_version(
        wildcards.sample_id,
        True,
        "0",
        lambda x: x >= 4,
    )

    return "xenium_ranger_new_cell_id" if is_xr_v4 else "xr_cell_id"


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
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        segmentation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0]
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadSegmentation2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        segmentation_id=r"(?!proseg|bats).*"
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**2 * 2048, 512000)
    script:
        "../../scripts/_standard_seurat_analysis/load_segmentation2seurat.R"

rule loadProseg2Seurat:
    input:
        unpack(get_input2_or_params4loadProseg2Seurat)
    output:
        protected(f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds')
    params:
        input_data=lambda wildcards: get_input2_or_params4loadProseg2Seurat(
            wildcards,
            for_input=False,
        ),
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        sample_id=lambda wildcards: wildcards.sample_id,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        control_gene_pat=sec.XENIUM_CONTROL_GENE_PAT,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        segmentation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0],
        xr_cell_id_col_name=get_xenium_cell_id4mapping,
        proseg_cell_id_col_name=get_imported_cell_id4mapping,
        use_mode_counts=use_mode_counts4loadProseg2Seurat,
        use_mapping=get_mapping_param4loadProseg2Seurat
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadProseg2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        segmentation_id=r"proseg_\w+"
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**2 * 4096, 1024000)
    script:
        "../../scripts/_standard_seurat_analysis/load_proseg2seurat.R"

rule loadBats2Seurat:
    input:
        dir=get_input2_or_params4loadSegmentation2Seurat,
        expected_counts=f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/raw_results/expected_counts.parquet'
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
        control_gene_pat=sec.XENIUM_CONTROL_GENE_PAT,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        segmentation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0]
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadBats2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        segmentation_id="bats"
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**2 * 4096, 1024000)
    script:
        "../../scripts/_standard_seurat_analysis/load_bats2seurat.R"
