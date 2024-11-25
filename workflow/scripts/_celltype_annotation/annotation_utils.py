# Functions for annotations .smk

def get_path_to_query_4annotation(wildcards) -> str:
    """
    Generate the default or specific path for the query used in annotation.
    """
    ret = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/std_seurat_objects/preprocessed_seurat.rds'

    annotation_mode = extract_layers_from_experiments(wildcards.annotation_id, [4])[0]
    if annotation_mode != "single_cell":
        ret = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/{annotation_mode}/seurat.rds'
    
    return ret


def get_path_to_reference_4referencebased_annotation(wildcards) -> str:
    """
    Generate the path for the reference object for reference-based annotation.
    """
    reference_path = lambda wildcards: get_dict_value(
        config,
        "experiments",
        "_cell_type_annotation",
        extract_layers_from_experiments(wildcards.sample_id, [0])[0],
        extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
        "path"
    )
    return reference_path


def get_output_folder_4referencebased_annotation(wildcards) -> str:
    """
    Generate the output folder path for storing reference-based annotation results.
    """
    output_folder = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/cell_type_annotation/{wildcards.annotation_id}'
    return output_folder
