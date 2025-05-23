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
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            "ovrlpy",
            f'signal_integrity_threshold={get_dict_value(
                config,
                "count_correction",
                "ovrlpy",
                "signal_integrity_threshold",
            )}',
            suffix,
        )
    elif re.match(
        r"^resolvi_unsupervised$",
        wildcards.count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            "resolvi_unsupervised",
            f'mixture_k={get_dict_value(
                config,
                "count_correction",
                "resolvi",
                "train",
                "mixture_k",
            )}',
            f'num_samples={get_dict_value(
                config,
                "count_correction",
                "resolvi",
                "predict",
                "num_samples",
            )}',
            suffix,
        )
    elif re.match(
        r"^resolvi_supervised$",
        wildcards.count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            wildcards.normalisation_id,
            wildcards.annotation_id,
            "resolvi_supervised",
            f'mixture_k={get_dict_value(
                config,
                "count_correction",
                "resolvi",
                "train",
                "mixture_k",
            )}',
            f'num_samples={get_dict_value(
                config,
                "count_correction",
                "resolvi",
                "predict",
                "num_samples",
            )}',
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

rule adaptCorrectedCountsBySplit2Seurat:
    input:
        raw_obj=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/raw_seurat.rds',
        corrected_counts=get_corrected_counts4adaptCorrectedCounts2Seurat
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/raw_seurat.rds')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/logs/adaptCorrectedCountsBySplit2Seurat.log'
    container:
        config["containers"]["r"]
    wildcard_constraints:
        count_correction_id=r"split.+"
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt**3 * 2048, 512000)
    script:
        "../../scripts/_standard_seurat_analysis/adapt_corrected_counts2seurat.R"

use rule adaptCorrectedCountsBySplit2Seurat as adaptCorrectedCountsByOvrlpy2Seurat with:
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/raw_seurat.rds')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/logs/adaptCorrectedCountsByOvrlpy2Seurat.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule adaptCorrectedCountsBySplit2Seurat as adaptCorrectedCountsByResolviUnsupervised2Seurat with:
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/raw_seurat.rds')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/logs/adaptCorrectedCountsByResolviUnsupervised2Seurat.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised"

use rule adaptCorrectedCountsBySplit2Seurat as adaptCorrectedCountsByResolviSupervised2Seurat with:
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/raw_seurat.rds')
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/logs/adaptCorrectedCountsByResolviSupervised2Seurat.log'
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised"
