import re


######################################
#              Subrules              #
######################################

include: '_segmentation/10x.smk'
include: '_segmentation/baysor.smk'


#######################################
#              Functions              #
#######################################

def get_input2getSegmentation(wildcards) -> str:
    if re.match(r"^10x_\d+?um$", wildcards.segmentation_id, flags=re.IGNORECASE):
        return f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/raw_results'
    elif wildcards.segmentation_id in ["baysor", "proseg"]:
        return f'{config["output_path"]}/segmentation/{wildcards.segmentation_id}/{wildcards.sample_id}/normalised_results'
    else:
        raise RuntimeError(f"Error! This workflow does not support method {wildcards.segmentation_id} for segmentation for now.")


#######################################
#                Rules                #
#######################################

rule getSegmentation:
    input:
        results_dir=lambda wildcards: get_input2getSegmentation(wildcards)
    output:
        directory(f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/results')
    log:
        f'{config["output_path"]}/segmentation/{{segmentation_id}}/{{sample_id}}/logs/getSegmentation.log'
    params:
        abs_input=lambda wildcards, input: os.path.abspath(
            input["results_dir"]
        )
    shell:
        "ln -s {params.abs_input} {output} && "
        "echo 'Done soft linking' > {log}"
