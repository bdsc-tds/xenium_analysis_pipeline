#######################################
#                Rules                #
#######################################

rule runResolviSupervisedTrain:
    input:
        unpack(get_seg_data4input2_or_param4runResolvi),
        annotation=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'
    output:
        directory(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/model')
    params:
        lambda wildcards: get_params4runResolvi(
            wildcards,
            for_training=True,
        )
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/logs/runResolviSupervisedTrain.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_mem_mb4runResolvi(
                wildcards,
                attempt,
                1,
            ),
            1024000,
        ),
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        slurm_extra=get_slurm_extra
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/resolvi_sample_training.py "
        "--path {params[0][data_dir]} "
        "--cell_type_labels {input.annotation} "
        "--out_dir_resolvi_model {output} "
        "--min_counts {params[0][min_counts]} "
        "--min_features {params[0][min_features]} "
        "--max_counts {params[0][max_counts]} "
        "--max_features {params[0][max_features]} "
        "--min_cells {params[0][min_cells]} "
        "--max_epochs {params[0][max_epochs]} "
        "--mixture_k {params[0][mixture_k]} "
        "-l {log}"

rule runResolviSupervisedPredict:
    input:
        unpack(get_seg_data4input2_or_param4runResolvi),
        model_dir=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/model',
        annotation=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'
    output:
        corrected_counts=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/corrected_counts.h5'),
        proportions=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/proportions.parquet')
    params:
        lambda wildcards: get_params4runResolvi(
            wildcards,
            for_training=False,
        )
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/logs/runResolviSupervisedPredict.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_mem_mb4runResolvi(
                wildcards,
                attempt**2,
                1,
            ),
            1024000
        ),
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        slurm_extra=get_slurm_extra
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/resolvi_sample_inference.py "
        "--path {params[0][data_dir]} "
        "--dir_resolvi_model {input.model_dir} "
        "--cell_type_labels {input.annotation} "
        "--out_file_resolvi_corrected_counts {output.corrected_counts} "
        "--out_file_resolvi_proportions {output.proportions} "
        "--min_counts {params[0][min_counts]} "
        "--min_features {params[0][min_features]} "
        "--max_counts {params[0][max_counts]} "
        "--max_features {params[0][max_features]} "
        "--min_cells {params[0][min_cells]} "
        "--num_samples {params[0][num_samples]} "
        "--batch_size {params[0][batch_size]} "
        "-l {log}"
