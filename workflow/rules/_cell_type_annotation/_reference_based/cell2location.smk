#######################################
#                Rules                #
#######################################

rule runCell2Location:
    input:
        query=get_path2query4annotation,
        reference=get_path2reference4reference_based_annotation
    output:
        h5ad=protected(f'{config["output_path"]}/cell2location/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.h5ad'),
        labels=protected(f'{config["output_path"]}/cell2location/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet'),
        scores=protected(f'{config["output_path"]}/cell2location/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet'),
        cell_abundance=protected(f'{config["output_path"]}/cell2location/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/cell_abundance.parquet'),
    params:
        annotation_id=lambda wildcards: wildcards.annotation_id,
        xe_assay=lambda wildcards: get_assay_name4annotation(wildcards, False),
        ref_assay=lambda wildcards: get_assay_name4annotation(wildcards, True),
        max_epochs_train=lambda wildcards: get_dict_value(
            config,
            "cell2location",
            "train",
            "max_epochs",
            replace_none=30000,
        ),
        max_epochs_predict=lambda wildcards: get_dict_value(
            config,
            "cell2location",
            "predict",
            "max_epochs",
            replace_none=1000,
        ),
        N_cells_per_location=lambda wildcards: get_dict_value(
            config,
            "cell2location",
            "predict",
            "N_cells_per_location",
            replace_none=5,
        ),
        detection_alpha=lambda wildcards: get_dict_value(
            config,
            "cell2location",
            "predict",
            "detection_alpha",
            replace_none=200,
        ),
    wildcard_constraints:
        annotation_id=r"reference_based/.+/cell2location/.+"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 10, 20480),
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        slurm_extra=get_slurm_extra
    log:
        f'{config["output_path"]}/cell2location/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runCell2Location.log'
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n cell2location python3 workflow/scripts/_cell_type_annotation/_reference_based/run_cell2location.py "
        "--query {input.query} "
        "--reference {input.reference} "
        "--annotation_id {params.annotation_id} "
        "--xe_assay {params.xe_assay} "
        "--ref_assay {params.ref_assay} "
        "--out_h5ad {output.h5ad} "
        "--out_labels {output.labels} "
        "--out_scores {output.scores} "
        "--out_cell_abundance {output.cell_abundance} "
        "--max_epochs_train {params.max_epochs_train} "
        "--max_epochs_predict {params.max_epochs_predict} "
        "--N_cells_per_location {params.N_cells_per_location} "
        "--detection_alpha {params.detection_alpha} "
        "-l {log}"
