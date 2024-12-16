#######################################
#                Rules                #
#######################################

rule runReferenceBasedRCTD:
    input:
      query=get_path2query4annotation,
      reference=get_path2reference4reference_based_annotation
    output:
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/output.rds'),
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/labels.csv'),
      protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/scores.csv')
    params:
      annotation_id=lambda wildcards: wildcards.annotation_id,
      #class_level = "Level2", #TODO: think how to provide class_level (Level{min(i-1, 1)}) for Level{i}. Extract levels hierarchi from experiments.yml?              
      ref_default_assay=cac.REF_SEURAT_DEFAULT_ASSAY,
      xe_default_assay=cac.XE_SEURAT_DEFAULT_ASSAY,
      REF_MIN_UMI=cac.REF_MIN_UMI,
      REF_MAX_UMI=cac.REF_MAX_UMI,
      XE_MIN_UMI=cac.XE_MIN_UMI,
      XE_MIN_counts=cac.XE_MIN_counts,
      CELL_MIN_INSTANCE=cac.CELL_MIN_INSTANCE, # RCTD-specific -> move to config or experiment
      UMI_min_sigma=100 # RCTD-specific -> move to config or experiment. Replace with one from `config -> ... -> rctd -> _mode -> _other_options` (shoould be increased for 5k panel)
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd/.+"
    threads:
        lambda wildcards: get_dict_value(
            config,
            "cell_type_annotation",
            "reference_based",
            "rctd",
            extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
            "_threads",
            replace_none=20
        )
    resources:
        mem_mb=lambda wildcards, input: max((os.path.getsize(input["query"]) >> 20) * 50, 20480)
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/cell_type_annotation/{{annotation_id}}/logs/runReferenceBasedRCTD.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_celltype_annotation/_reference_based/rctd.R"
