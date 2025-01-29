#######################################
#              Functions              #
#######################################

def get_input2_or_params4loadSegmentation2Seurat(wildcards, for_input: bool = True) -> str:
    ret = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/normalised_results'

    if for_input:
        return ret

    return normalise_path(
        ret,
        candidate_paths=("outs",),
        pat_flags=re.IGNORECASE,
        return_dir=True,
        check_exist=False
    )


#######################################
#                Rules                #
#######################################

rule loadSegmentation2Seurat:
    input:
        get_input2_or_params4loadSegmentation2Seurat
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/std_seurat_objects/raw_seurat.rds')
    params:
        data_dir=lambda wildcards: get_input2_or_params4loadSegmentation2Seurat(
            wildcards,
            for_input=False
        ),
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        sample_id=lambda wildcards: wildcards.sample_id,
        segmentation_id=lambda wildcards: wildcards.segmentation_id
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/loadSegmentation2Seurat.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/load_segmentation2seurat.R"
