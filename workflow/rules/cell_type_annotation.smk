import scripts._cell_type_annotation.cell_type_annotation_constants as cac


#######################################
#              Functions              #
#######################################

def get_path2query4annotation(wildcards, is_post_correction: bool = False) -> str:
    """
    Get the default or specific path to the query used in annotation.
    """
    suffix: str = "preprocessed/preprocessed_seurat.rds"

    annotation_mode = extract_layers_from_experiments(wildcards.annotation_id, [4])[0]

    if annotation_mode == "single_cell":
        if not is_post_correction:
            return os.path.join(
                f'{config["output_path"]}/std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}',
                suffix,
            )
        
        _prefix: str = f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}'

        if re.fullmatch(
            r"^ovrlpy$",
            wildcards.count_correction_id,
            flags=re.IGNORECASE,
        ) is not None:
            return os.path.join(
                _prefix,
                "ovrlpy",
                f'signal_integrity_threshold={get_dict_value(
                    config,
                    "count_correction",
                    "ovrlpy",
                    "signal_integrity_threshold",
                )}',
                wildcards.normalisation_id,
                suffix,
            )
        elif re.fullmatch(
            r"^resolvi_unsupervised$",
            wildcards.count_correction_id,
            flags=re.IGNORECASE,
        ) is not None:
            return os.path.join(
                _prefix,
                "resolvi_unsupervised",
                f'mixture_k={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "train",
                    "mixture_k",
                )}',
                f'num_samples={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "predict",
                    "num_samples",
                )}',
                wildcards.normalisation_id,
                suffix,
            )
        elif re.fullmatch(
            r"^resolvi_supervised$",
            wildcards.count_correction_id,
            flags=re.IGNORECASE,
        ) is not None:
            return os.path.join(
                _prefix,
                wildcards.normalisation_id,
                wildcards.annotation_id,
                "resolvi_supervised",
                f'mixture_k={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "train",
                    "mixture_k",
                )}',
                f'num_samples={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "predict",
                    "num_samples",
                )}',
                wildcards.normalisation_id,
                suffix,
            )
        else:
            return os.path.join(
                _prefix,
                wildcards.normalisation_id,
                wildcards.annotation_id,
                wildcards.count_correction_id,
                wildcards.normalisation_id,
                suffix,
            )
    else:
        raise RuntimeError(f"Error! Unknown mode for annotation: {annotation_mode}. Valid modes include 'single_cell'.")


def get_path2reference4reference_based_annotation(wildcards) -> str:
    """
    Get the path to the reference object for reference-based annotation.
    """
    return get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_REFERENCE_FILES_NAME,
        extract_layers_from_experiments(wildcards.sample_id, [0])[0],
        extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
    )


def get_assay_name4annotation(wildcards, for_ref: bool) -> str:
    if wildcards.normalisation_id == "lognorm":
        return cac.REF_SEURAT_DEFAULT_ASSAY if for_ref else cac.XE_SEURAT_DEFAULT_ASSAY
    elif wildcards.normalisation_id == "sctransform":
        return cac.REF_SEURAT_ALT_ASSAY if for_ref else cac.XE_SEURAT_ALT_ASSAY

    raise RuntimeError(f"Error! Unknown normalisation method: {wildcards.normalisation_id}. Valid methods include 'lognorm' and 'sctransform'.")


# Shared by rctd.smk (R) and rctd_py.smk (python)
def get_class_level4runReferenceBasedRCTD(wildcards) -> str | None:
    method_name: str = extract_layers_from_experiments(wildcards.annotation_id, [2])[0]

    if method_name == "rctd_class_unaware":
        return None
    elif method_name == "rctd_class_aware":
        current_level: str = extract_layers_from_experiments(wildcards.annotation_id, [3])[0]
        available_levels: list[str] = get_dict_value(
                config,
                "experiments",
                cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
                cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_LEVELS_NAME,
                extract_layers_from_experiments(wildcards.sample_id, [0])[0],
                extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
            )

        try:
            current_level_idx: int = available_levels.index(current_level)
        except ValueError:
            raise RuntimeError(f"Error! Cannot find {current_level} in [{",".join(available_levels)}].")
        
        return available_levels[0] if current_level_idx == 0 else available_levels[current_level_idx - 1]
    else:
        raise RuntimeError(f"Error! Unknown method name: {method_name}.")


######################################
#              Subrules              #
######################################

# RCTD implementation from cell_type_annotation.reference_based.rctd.language: R (spacexr, default) or
# python (rctd-py); only the selected one is included, both write the same output files
_RCTD_LANGUAGE: set[str] = {
    str(v.get("language", "R"))
    for v in (get_dict_value(config, "cell_type_annotation", "reference_based", "rctd", inexist_key_ok=True) or {}).values()
    if isinstance(v, dict)
}

if _RCTD_LANGUAGE == {"python"}:
    include: "_cell_type_annotation/_reference_based/rctd_py.smk"
elif _RCTD_LANGUAGE <= {"R"}:
    include: "_cell_type_annotation/_reference_based/rctd.smk"
else:
    raise RuntimeError(f"Error! cell_type_annotation.reference_based.rctd.language must be 'R' or 'python', got {_RCTD_LANGUAGE}.")

include: "_cell_type_annotation/_reference_based/singler.smk"
include: "_cell_type_annotation/_reference_based/seurat.smk"
include: "_cell_type_annotation/_reference_based/xgboost.smk"
