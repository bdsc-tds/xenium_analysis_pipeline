#######################################
#              Functions              #
#######################################

def get_seg_data4input2_or_param4runResolvi(wildcards, for_input: bool = True) -> str | list[str]:
    prefix: str = f'{config["output_path"]}/segmentation'
    ret: str = ""

    if re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            f'proseg/{wildcards.sample_id}/raw_results'
        )

        if for_input:
            ret = [
                os.path.join(
                    ret,
                    i
                ) for i in [
                    "cell-metadata.csv.gz",
                    "transcript-metadata.csv.gz",
                    "cell-polygons.geojson.gz",
                    "expected-counts.csv.gz"
                ]
            ]
    else:
        ret = os.path.join(
            prefix,
            f"{wildcards.segmentation_id}/{wildcards.sample_id}/normalised_results",
        )

        if not for_input:
            ret = normalise_path(
                ret,
                candidate_paths=("outs",),
                pat_anchor_file=r"transcripts\.parquet",
                pat_flags=re.IGNORECASE,
                return_dir=True,
                check_exist=False
            )

    return ret


def get_params4runResolvi(wildcards, for_training: bool) -> dict[str, Any]:
    ret: dict[atr, Any] = dict()

    ret["data_dir"] = get_seg_data4input2_or_param4runResolvi(
        wildcards,
        for_input=False,
    )

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
    ret: int = get_size(
        get_seg_data4input2_or_param4runResolvi(
            wildcards,
        ),
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

include: 'resolvi_unsupervised.smk'
include: 'resolvi_supervised.smk'
