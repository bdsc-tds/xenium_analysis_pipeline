#######################################
#              Functions              #
#######################################

def get_corrected_counts4adaptCorrectedCounts2Seurat(wildcards) -> str:
    prefix: str = f'{config["output_path"]}/count_correction/{wildcards.segmentation_id}/{wildcards.sample_id}'
    suffix: str = "corrected_counts.h5"

    ret: str = ""

    if re.match(
        r"^ovrlpy$",
        wildcards.count_correction_id,
        flags=re.IGNOREXASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            "ovrlpy",
            f'signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}',
            suffix,
        )
    elif re.match(
        r"^resolvi_unsupervised$",
        wildcards.count_correction_id,
        flags=re.IGNOREXASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            "resolvi_unsupervised",
            suffix,
        )
    else:
        ret = os.path.join(
            prefix,
            wildcards.normalisation_id,
            wildcards.annotation_id,
            wildcards.count_correction_id,
            suffix,
        )

    return ret


#######################################
#                Rules                #
#######################################

rule adaptCorrectedCounts2Seurat:
    input:
        raw_obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds',
        corrected_counts=get_corrected_counts4adaptCorrectedCounts2Seurat
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/raw_seurat.rds')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/logs/adaptCorrectedCounts2Seurat.log'
    container:
        config["containers"]["r"]
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**2 * 2048, 512000)
    script:
        "../../scripts/_standard_seurat_analysis/adapt_corrected_counts2seurat.R"
