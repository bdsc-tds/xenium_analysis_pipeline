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
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/normalised_results'
    output:
        protected(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/normalised_results/_auxiliary_files.tar')
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/zipNormalisedAuxiliary10xFiles.log'
    params:
        abs_output=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/normalised_results/_auxiliary_files.tar'
        ),
        filenames=get_auxiliary_10x_files
    shell:
        "cd {input} && "
        "tar --remove-files -cf {params.abs_output} {params.filenames} &> {log}"
