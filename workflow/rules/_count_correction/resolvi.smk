#######################################
#              Functions              #
#######################################

def get_seg_data4input2_or_param4runResolvi(wildcards, for_input: bool = True) -> dict[str, str]:
    prefix: str = f'{config["output_path"]}/segmentation'
    ret: dict[str, str] = {}

    if re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        dir_path: str = os.path.join(
            prefix,
            f'proseg/{wildcards.sample_id}/raw_results'
        )

        if for_input:
            ret["cell_metadata"] = os.path.join(
                dir_path,
                "cell-metadata.csv.gz",
            )
            ret["transcript_metadata"] = os.path.join(
                dir_path,
                "transcript-metadata.csv.gz",
            )
            ret["cell_polygons"] = os.path.join(
                dir_path,
                "cell-polygons.geojson.gz",
            )
            ret["expected_counts"] = os.path.join(
                dir_path,
                "expected-counts.csv.gz",
            )
        else:
            ret["data_dir"] = dir_path
    else:
        seg_id: str = "proseg" if re.match(
            r"^proseg_mode$",
            wildcards.segmentation_id,
            flags=re.IGNORECASE,
        ) is not None else wildcards.segmentation_id

        dir_path = os.path.join(
            prefix,
            f"{seg_id}/{wildcards.sample_id}/normalised_results",
        )

        ret["data_dir"] = dir_path if for_input else normalise_path(
            dir_path,
            candidate_paths=("outs",),
            pat_anchor_file=r"transcripts.parquet",
            pat_flags=re.IGNORECASE,
            return_dir=True,
            check_exist=False
        )

    return ret


def get_params4runResolvi(wildcards, for_training: bool) -> dict[str, Any]:
    ret: dict[atr, Any] = dict()

    ret = get_seg_data4input2_or_param4runResolvi(
        wildcards,
        for_input=False,
    )
    assert "data_dir" in ret

    ret["min_counts"] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        extract_layers_from_experiments(
            wildcards.sample_id,
            [0, 1],
        )[0],
        "min_counts",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_counts",
            replace_none=10,
        ),
        inexist_key_ok=True,
    )

    ret["min_features"] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        extract_layers_from_experiments(
            wildcards.sample_id,
            [0, 1],
        )[0],
        "min_features",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_features",
            replace_none=5,
        ),
        inexist_key_ok=True,
    )

    ret["max_counts"] = float(get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        extract_layers_from_experiments(
            wildcards.sample_id,
            [0, 1],
        )[0],
        "max_counts",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "max_counts",
            replace_none="Inf",
        ),
        inexist_key_ok=True,
    ))

    ret["max_features"] = float(get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        extract_layers_from_experiments(
            wildcards.sample_id,
            [0, 1],
        )[0],
        "max_features",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "max_features",
            replace_none="Inf",
        ),
        inexist_key_ok=True,
    ))

    ret["min_cells"] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        extract_layers_from_experiments(
            wildcards.sample_id,
            [0, 1],
        )[0],
        "min_cells",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_cells",
            replace_none=5,
        ),
        inexist_key_ok=True,
    )

    if for_training:
        ret["max_epochs"] = get_dict_value(
            config,
            "count_correction",
            wildcards.count_correction_id,
            "train",
            "max_epochs",
            replace_none=50,
        )

        ret["mixture_k"] = get_dict_value(
            config,
            "count_correction",
            wildcards.count_correction_id,
            "train",
            "mixture_k",
            replace_none=50,
        )
    else:
        ret['num_samples'] = get_dict_value(
            config,
            "count_correction",
            wildcards.count_correction_id,
            "predict",
            "num_samples",
            replace_none=30,
        )

        ret["batch_size"] = get_dict_value(
            config,
            "count_correction",
            wildcards.count_correction_id,
            "predict",
            "batch_size",
            replace_none=1000,
        )

    return ret


def get_mem_mb4runResolvi(wildcards, attempt, multiplier: int = 1) -> int:
    ret: int = sum(
        [
            get_size(i)
            for i in get_seg_data4input2_or_param4runResolvi(
                wildcards,
                for_input=False,
            ).values()
        ]
    ) * 1e-6 * attempt

    if re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret *= 10

    if not _use_gpu():
        ret *= 2

    ret *= multiplier
    return ret


######################################
#              Subrules              #
######################################

include: '_resolvi/resolvi_unsupervised.smk'
include: '_resolvi/resolvi_supervised.smk'
