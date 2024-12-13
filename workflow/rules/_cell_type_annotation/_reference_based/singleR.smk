#######################################
#              Functions              #
#######################################

#include: "../../../scripts/_celltype_annotation/annotation_utils.py" # @Senbai, any better way? -> moved to parent .smk (~/workflow/rules/celltype_annotation.smk)


#######################################
#                Rules                #
#######################################


rule runSingleR:
    input:
      query = get_path_to_query_4annotation,
      reference = get_path_to_reference_4referencebased_annotation, 
    output:
      protected(f'{get_output_folder_4referencebased_annotation}/output.rds'),
      protected(f'{get_output_folder_4referencebased_annotation}/labels.csv'),
      protected(f'{get_output_folder_4referencebased_annotation}/scores.csv')
    params:
      annotation_id = lambda wildcards: wildcards.annotation_id,
      ref_default_assay = sac.REF_SEURAT_DEFAULT_ASSAY, 
      xe_default_assay = sac.XE_SEURAT_DEFAULT_ASSAY, 
      REF_MIN_UMI = sac.REF_MIN_UMI, 
      REF_MAX_UMI = sac.REF_MAX_UMI, 
      XE_MIN_UMI = sac.XE_MIN_UMI,
      XE_MIN_counts = sac.XE_MIN_counts,
      # singleR-specific params
      singleR_genes = sac.singleR_genes,
      de_method = sac.singleR_de_method, 
      aggr_ref = sac.singleR_aggr_ref, 
      aggr_args = sac.singleR_aggr_args # list(rank = 50, power = 0.7) # not sure how to 
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/{{annotation_id}.replace("/", "_")}.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_celltype_annotation/_reference_based/singleR.R"
