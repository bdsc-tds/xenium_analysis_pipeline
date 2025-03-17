import scripts._cell_type_annotation.cell_type_annotation_constants as cac


#######################################
#              Functions              #
#######################################

def get_path2query4annotation(wildcards, is_post_correction: bool = False) -> str:
    """
    Get the default or specific path to the query used in annotation.
    """
    annotation_mode = extract_layers_from_experiments(wildcards.annotation_id, [4])[0]

    if annotation_mode == "single_cell":
        if is_post_correction:
            return f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}/{wildcards.annotation_id}/{wildcards.count_correction_id}/{wildcards.normalisation_id}/preprocessed/preprocessed_seurat.rds'
        else:
            return f'{config["output_path"]}/std_seurat_analysis/{wildcards.segmentation_id}/{wildcards.sample_id}/{wildcards.normalisation_id}/preprocessed/preprocessed_seurat.rds'
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


######################################
#              Subrules              #
######################################

include: "_cell_type_annotation/_reference_based/rctd.smk"
include: "_cell_type_annotation/_reference_based/singler.smk"
include: "_cell_type_annotation/_reference_based/seurat.smk"
include: "_cell_type_annotation/_reference_based/xgboost.smk"
