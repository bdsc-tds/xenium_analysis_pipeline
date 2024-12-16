import scripts._cell_type_annotation.cell_type_annotation_constants as cac


#######################################
#              Functions              #
#######################################

def get_path2query4annotation(wildcards) -> str:
    """
    Get the default or specific path to the query used in annotation.
    """
    annotation_mode = extract_layers_from_experiments(wildcards.annotation_id, [4])[0]

    if annotation_mode == "single_cell":
        return f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/std_seurat_objects/preprocessed_seurat.rds'
    elif annotation_mode == "cluster": # TODO
        return f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/{annotation_mode}/seurat.rds'
    else:
        raise RuntimeError(f"Error! Unknown mode for annotation: {annotation_mode}. Valid modes include 'single_cell' and 'cluster'.")


def get_path2reference4reference_based_annotation(wildcards) -> str:
    """
    Get the path to the reference object for reference-based annotation.
    """
    return get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_REFERENCE_FILE_NAME,
        extract_layers_from_experiments(wildcards.sample_id, [0])[0],
        extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0]
    )


######################################
#              Subrules              #
######################################

include: "_cell_type_annotation/_reference_based/rctd.smk"
include: "_cell_type_annotation/_reference_based/singler.smk"
# include: "_cell_type_annotation/_reference_based/xgboost.smk"
# include: "_cell_type_annotation/_reference_based/seurat_transfer.smk"
# include: "_cell_type_annotation/_reference_based/tangram.smk" # placeholder
