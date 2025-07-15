#######################################
#              Functions              #
#######################################

# level: 0 for condition-wise, 1 for gene-panel-wise, 2 for donor-wise
def get_input2_or_params4gatherSegmentedSamples(
    wildcards,
    *,
    level: int,
    for_input: bool = True,
) -> list[str] | dict[str, str]:
    if level == 0:
        level_name: str = cc.EXPERIMENTS_COLLECTIONS_CONDITIONS_NAME
        level_id: str = wildcards.condition_id
    elif level == 1:
        level_name = cc.EXPERIMENTS_COLLECTIONS_GENE_PANELS_NAME
        level_id = wildcards.gene_panel_id
    elif level == 2:
        level_name = cc.EXPERIMENTS_COLLECTIONS_DONORS_NAME
        level_id = wildcards.donor_id
    else:
        raise ValueError("Error! Level must be 0, 1, or 2")

    prefix: str = f'{config["output_path"]}/joint_scanpy_analysis'

    samples: list[str] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_COLLECTIONS_NAME,
        level_name,
        level_id,
    )

    ret: dict[str: str] = {
        "_".join(
            [
                wildcards.segmentation_id,
                i,
            ]
        ): os.path.abspath(
            os.path.join(
                prefix,
                wildcards.segmentation_id,
                i,
                "raw_spatialdata.h5ad",
            )
        ) for i in samples
    }

    if for_input:
        return list(ret.values())
    else:
        return ret


# level: 0 for condition-wise, 1 for gene-panel-wise, 2 for donor-wise
def get_params4runAnalysis(
    wildcards,
    level: int,
) -> dict[str, Any]:
    if level == 0:
        level_name: str = "condition_wise"
    elif level == 1:
        level_name = "gene_panel_wise"

        gene_panel_id: str = extract_layers_from_experiments(
            wildcards.gene_panel_id,
            [0, 1],
        )[0]
    elif level == 2:
        level_name = "donor_wise"
    else:
        raise ValueError("Error! Level must be 0, 1, or 2")

    ret: dict[str, Any] = {}

    ret["use_gpu"] = "--use_gpu" if _use_gpu() else ""
    ret["n_comps"] = get_dict_value(
        config,
        "joint_scanpy_analysis",
        level_name,
        "n_comps",
        replace_none=50,
    )

    ret["n_neighbors"] = get_dict_value(
        config,
        "joint_scanpy_analysis",
        level_name,
        "n_neighbors",
        replace_none=50,
    )

    ret["min_dist"] = get_dict_value(
        config,
        "joint_scanpy_analysis",
        level_name,
        "min_dist",
        replace_none=0.3,
    )

    ret["metric"] = get_dict_value(
        config,
        "joint_scanpy_analysis",
        level_name,
        "metric",
        replace_none="cosine",
    )

    ret["min_counts"]=get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        gene_panel_id,
        "min_counts",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_counts",
            replace_none=10,
        ),
        inexist_key_ok=True,
    ) if level == 1 else get_dict_value(
        config,
        "standard_seurat_analysis",
        "qc",
        "min_counts",
        replace_none=10,
    )

    ret["min_features"] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        gene_panel_id,
        "min_features",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_features",
            replace_none=5,
        ),
        inexist_key_ok=True,
    ) if level == 1 else get_dict_value(
        config,
        "standard_seurat_analysis",
        "qc",
        "min_features",
        replace_none=5,
    )

    ret["max_counts"] = float(get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        gene_panel_id,
        "max_counts",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "max_counts",
            replace_none="Inf",
        ),
        inexist_key_ok=True,
    ) if level == 1 else get_dict_value(
        config,
        "standard_seurat_analysis",
        "qc",
        "max_counts",
        replace_none="Inf",
    ))
    
    ret["max_features"] = float(get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        gene_panel_id,
        "max_features",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "max_features",
            replace_none="Inf",
        ),
        inexist_key_ok=True,
    ) if level == 1 else get_dict_value(
        config,
        "standard_seurat_analysis",
        "qc",
        "max_features",
        replace_none="Inf",
    ))

    ret["min_cells"] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
        gene_panel_id,
        "min_cells",
        replace_none=get_dict_value(
            config,
            "standard_seurat_analysis",
            "qc",
            "min_cells",
            replace_none=5,
        ),
        inexist_key_ok=True,
    ) if level == 1 else get_dict_value(
        config,
        "standard_seurat_analysis",
        "qc",
        "min_cells",
        replace_none=5,
    )

    return ret


# level: 0 for condition-wise, 1 for gene-panel-wise, 2 for donor-wise
def get_mem_mb4runAnalysis(
    wildcards,
    attempt,
    level: int,
    multiplier: int = 1
) -> int:
    ret: int = sum(
        get_size(i)
        for i in get_input2_or_params4gatherSegmentedSamples(
            wildcards,
            level=level,
        )
    ) * 1e-6 * attempt
    ret *= multiplier
    return ret


######################################
#              Subrules              #
######################################

include: "_joint_scanpy_analysis/load_segmentation2spatialdata.smk"
include: "_joint_scanpy_analysis/gene_panel_wise_analysis.smk"
include: "_joint_scanpy_analysis/condition_wise_analysis.smk"
