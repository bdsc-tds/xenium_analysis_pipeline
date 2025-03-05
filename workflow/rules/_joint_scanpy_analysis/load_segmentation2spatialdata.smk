#######################################
#              Functions              #
#######################################

def get_input2loadSegmentation2SpatialData(wildcards) -> dict[str, str]:
    ret: dict[str, str] = get_seg_data4input2_or_param4runResolvi(
        wildcards,
        for_input=True,
    )

    if re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret["mapping"] = f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/mapped_cell_ids/mapped_cell_ids.parquet'

    return ret


def get_cmd_args4loadSegmentation2SpatialData(wildcards, input) -> dict[str, str]:
    data_dir: dict[str, str] = get_seg_data4input2_or_param4runResolvi(
        wildcards,
        for_input=False,
    )
    assert "data_dir" in data_dir

    args: list[str] = [
        "-i",
        data_dir["data_dir"],
        "--sample_id",
        wildcards.sample_id,
        "--segmentation_id",
        wildcards.segmentation_id,
        "--condition",
        extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        "--gene_panel",
        extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        "--donor",
        extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        "--sample",
        extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        "--segmentation_method",
        extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0],
    ]

    if re.match(
        r"^proseg_\w+",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        assert "mapping" in input
        args.extend(
            [
                "--in_mapping",
                input["mapping"],
                "--cell_id_col_name",
                "xr_cell_id" if use_mode_counts4loadProseg2Seurat(
                    wildcards,
                ) else "proseg_cell_id",
            ]
        )

        if get_mapping_param4loadProseg2Seurat(
            wildcards,
        ):
            args.append(
                "--use_mapping",
            )

    return " ".join(args)


######################################
#                Rules                #
#######################################

rule loadSegmentation2SpatialData:
    input:
        unpack(get_input2loadSegmentation2SpatialData)
    output:
        protected(f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{sample_id}}/raw_spatialdata.h5ad')
    params:
        cmd_args=lambda wildcards, input: get_cmd_args4loadSegmentation2SpatialData(
            wildcards,
            input,
        )
    log:
        f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{sample_id}}/logs/loadSegmentation2SpatialData.log'
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: max(
            get_size(
                get_seg_data4input2_or_param4runResolvi(
                    wildcards,
                    for_input=False,
                )["data_dir"]
            ) * 10**-6 * attempt * 4,
            2048
        )
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_joint_scanpy_analysis/load_segmentation2spatialdata.py "
        "{params.cmd} "
        "-l {log}"
