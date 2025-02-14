#######################################
#              Functions              #
#######################################

def _use_gpu() -> bool:
    return get_dict_value(
        config,
        "gpu",
        "available",
        replace_none=False,
    )

def get_slurm_gpu_partition_name(wildcards) -> str:
    return get_dict_value(
        config,
        "gpu",
        "partition_name",
        replace_none="gpu",
    )

def get_slurm_extra(wildcards, num_device: int = 1) -> str:
    return f"--gres=gpu:{num_device}" if _use_gpu() else ""


######################################
#              Subrules              #
######################################

include: '_segmentation/10x.smk'
include: '_segmentation/baysor.smk'
include: '_segmentation/proseg.smk'
include: '_segmentation/segger.smk'


#######################################
#                Rules                #
#######################################

rule zipNormalisedAuxiliary10xFiles:
    input:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/normalised_results'
    output:
        protected(f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/normalised_results/_auxiliary_files.tar')
    log:
        f'{config["output_path"]}/segmentation/{{compact_segmentation_id}}/{{sample_id}}/logs/zipNormalisedAuxiliary10xFiles.log'
    params:
        abs_output=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/{wildcards.compact_segmentation_id}/{wildcards.sample_id}/normalised_results/_auxiliary_files.tar'
        ),
        filenames=get_auxiliary_10x_files
    shell:
        "cd {input} && "
        "tar --remove-files -cf {params.abs_output} {params.filenames} &> {log}"
