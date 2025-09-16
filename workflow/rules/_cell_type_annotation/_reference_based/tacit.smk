#######################################
#                Rules                #
#######################################

rule runReferenceBasedTACIT:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
        rctd=path_2_RCTD_results # @Senbai, please get path to class-aware RCTD at same level and reference, eg `/users/mbilous/spatial/data/xenium_paper/xenium/processed/cell_type_annotation/10x_5um/NSCLC/lung/0PSV/0PSV/lognorm/reference_based/matched_reference_combo/rctd_class_aware/Level2.1/single_cell/output.rds`
    output:
        tacit_output=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds'),
        labels=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        tacit_weights=protected(f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/assignment_matrix.parquet')
    params:
        mode="binary", #"binary" or "continuous" # Need Senbai's help ><
        N_binary_markers_per_cell_type = 10,
        N_PCs = 50,
        r = 10, # minicluster size in TACIT
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        annotation_id=lambda wildcards: wildcards.annotation_id,
        ref_assay=lambda wildcards: get_assay_name4annotation(wildcards, True),
        xe_assay=lambda wildcards: get_assay_name4annotation(wildcards, False),
        REF_MIN_UMI=cac.REF_MIN_UMI,
        REF_MAX_UMI=cac.REF_MAX_UMI,
        XE_MIN_UMI=cac.XE_MIN_UMI,
        XE_MIN_counts=cac.XE_MIN_counts,
        cell_min_instance=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_CELL_MIN_INSTANCES_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0])[0],
            extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
            replace_none=25,
        )
    wildcard_constraints:
        annotation_id=r"reference_based/.+/tacit_.+"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 20, 20480)
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedTACIT.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/tacit.R"

