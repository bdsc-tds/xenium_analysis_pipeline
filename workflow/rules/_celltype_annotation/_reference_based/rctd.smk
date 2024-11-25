#######################################
#              Functions              #
#######################################

get_path_to_query_4annotation(wildcards) -> str: # Not RCTD-specific, move to parent folder or utils file
  # default path (to processed single-cell object)
  ret = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/std_seurat_objects/preprocessed_seurat.rds'
  
  annotation_mode = extract_layers_from_experiments(wildcards.annotation_id, [4])[0]
  if annotation_mode != "single_cell":
    ret = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/{annotation_mode}/seurat.rds' # yet to be generated
  
  return ret
      

get_path_to_reference_4referencebased_annotation(wildcards) -> str:
  
  reference_path = lambda wildcards: get_dict_value(
    config,
    "experiments",
    "_cell_type_annotation", # TODO: make constant @Senbai?
    extract_layers_from_experiments(wildcards.sample_id, [0])[0], # get disease 
    extract_layers_from_experiments(wildcards.annotation_id, [0,1])[0], # get annotation_approach == reference_based and reference_type # NOT SURE this function can be used with wildcards.annotation_id
    "path"
  )

  return reference_path
    

get_output_folder_4referencebased_annotation(wildcards) -> str:
  
  output_folder = f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/cell_type_annotation/{wildcards.annotation_id}'
  
  return output_folder


#######################################
#                Rules                #
#######################################


rule runRCTD:
    input:
      query = get_path_to_query_4annotation,
      reference = get_path_to_reference_4referencebased_annotation, 
    output:
      protected(f'{get_output_folder_4referencebased_annotation}/output.rds'),
      protected(f'{get_output_folder_4referencebased_annotation}/labels.csv'),
      protected(f'{get_output_folder_4referencebased_annotation}/scores.csv')
    params:
      annotation_id = lambda wildcards: wildcards.annotation_id,
      #class_level = "Level2", #TODO: think how to provide class_level (Level{min(i-1, 1)}) for Level{i}. Extract levels hierarchi from experiments.yml?              
      ref_default_assay = sac.REF_SEURAT_DEFAULT_ASSAY, 
      xe_default_assay = sac.XE_SEURAT_DEFAULT_ASSAY, 
      REF_MIN_UMI = sac.REF_MIN_UMI, 
      REF_MAX_UMI = sac.REF_MAX_UMI, 
      XE_MIN_UMI = sac.XE_MIN_UMI,
      XE_MIN_counts = sac.XE_MIN_counts,
      CELL_MIN_INSTANCE = sac.CELL_MIN_INSTANCE, # RCTD-specific -> move to config or experiment
      UMI_min_sigma = 100 # RCTD-specific -> move to config or experiment. Replace with one from `config -> ... -> rctd -> _mode -> _other_options`
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/{{annotation_id}.replace("/", "_")}.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_celltype_annotation/_reference_based/rctd.R"
