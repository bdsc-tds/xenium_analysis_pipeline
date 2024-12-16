#######################################
#                Rules                #
#######################################

rule runReferenceBasedSeurat:
    input:
      query=get_path2query4annotation,
      reference=get_path2reference4reference_based_annotation
    output:
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/output.rds'),
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/labels.csv'),
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/scores.csv')
    params:
      annotation_id=lambda wildcards: wildcards.annotation_id,
      ref_default_assay=cac.REF_SEURAT_DEFAULT_ASSAY,
      xe_default_assay=cac.XE_SEURAT_DEFAULT_ASSAY,
      REF_MIN_UMI=cac.REF_MIN_UMI,
      REF_MAX_UMI=cac.REF_MAX_UMI,
      XE_MIN_UMI=cac.XE_MIN_UMI,
      XE_MIN_counts=cac.XE_MIN_counts,
      dims=cac.dims # 1:50
    wildcard_constraints:
        annotation_id=r"reference_based/.+/seurat/.+"
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/logs/runReferenceBasedSeurat.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_celltype_annotation/_reference_based/seurat.R"
