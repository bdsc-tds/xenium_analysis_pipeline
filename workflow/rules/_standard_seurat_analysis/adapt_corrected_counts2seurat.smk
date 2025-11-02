#######################################
#              Functions              #
#######################################

def get_corrected_counts4adaptCorrectedCounts2Seurat(wildcards) -> str:
    suffix: str = "corrected_counts.h5"

    if re.match(
        r"^ovrlpy$",
        wildcards.count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        return os.path.join(
            config["output_path"],
            "count_correction",
            wildcards.segmentation_id4ovrlpy,
            wildcards.sample_id,
            "ovrlpy",
            f'signal_integrity_threshold={get_dict_value(
                config,
                "count_correction",
                "ovrlpy",
                "signal_integrity_threshold",
            )}',
            suffix,
        )

    _prefix: str = f'{config["output_path"]}/count_correction/{wildcards.segmentation_id}/{wildcards.sample_id}'

    if re.match(
        r"^resolvi_unsupervised$",
        wildcards.count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        return os.path.join(
            _prefix,
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
        return os.path.join(
            _prefix,
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
        return os.path.join(
            _prefix,
            wildcards.normalisation_id,
            wildcards.annotation_id,
            wildcards.count_correction_id,
            suffix,
        )

def get_inpu2adaptCorrectedCountsByOvrlpy2Seurat(wildcards) -> dict[str, str]:
    ret: dict[str, str] = {
        "corrected_counts": get_corrected_counts4adaptCorrectedCounts2Seurat(
            wildcards,
        ),
    }

    if re.fullmatch(
        r"^raw$",
        wildcards.segmentation_id4ovrlpy,
        flags=re.IGNORECASE,
    ) is not None:
        ret["raw_obj"] = f'{config["output_path"]}/std_seurat_analysis/{wildcards.segmentation_id4ovrlpy}/{wildcards.sample_id}/raw_seurat.rds'
    elif re.fullmatch(
        r"^proseg$",
        wildcards.segmentation_id4ovrlpy,
        flags=re.IGNORECASE,
    ) is not None:
        avail_proseg_id: list[str] = sorted(
            [
                i
                for i in SEGMENTATION_ID
                if re.fullmatch(
                    r"^proseg_.+$",
                    i,
                    flags=re.IGNORECASE,
                ) is not None
            ]
        )
        assert len(avail_proseg_id) > 0, "Error! No available proseg segmentation id found."

        ret["raw_obj"] = f'{config["output_path"]}/std_seurat_analysis/{avail_proseg_id[0]}/{wildcards.sample_id}/raw_seurat.rds'
    else:
        raise RuntimeError(f'Error! Unknown segmentation id: {wildcards.segmentation_id4ovrlpy}.')
    
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
    params:
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        segmentation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0],
        sample_id=lambda wildcards: wildcards.sample_id,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        normalisation_id=lambda wildcards: wildcards.normalisation_id,
        annotation_id=lambda wildcards: wildcards.annotation_id,
        annotation_approach=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            0,
        )[0],
        annotation_reference_name=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            1,
        )[0],
        annotation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            2,
        )[0],
        annotation_level=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            3,
        )[0],
        annotation_mode=lambda wildcards: extract_layers_from_experiments(
            wildcards.annotation_id,
            4,
        )[0],
        count_correction_id=lambda wildcards: wildcards.count_correction_id
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
    input:
        unpack(get_inpu2adaptCorrectedCountsByOvrlpy2Seurat)
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id4ovrlpy}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/raw_seurat.rds')
    params:
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        segmentation_id=lambda wildcards: wildcards.segmentation_id4ovrlpy,
        segmentation_method=lambda wildcards: wildcards.segmentation_id4ovrlpy,
        sample_id=lambda wildcards: wildcards.sample_id,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        normalisation_id=None,
        annotation_id=None,
        annotation_approach=None,
        annotation_reference_name=None,
        annotation_method=None,
        annotation_level=None,
        annotation_mode=None,
        count_correction_id=lambda wildcards: wildcards.count_correction_id
    log:
        f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id4ovrlpy}}/{{sample_id}}/{{count_correction_id}}/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/logs/adaptCorrectedCountsByOvrlpy2Seurat.log'
    wildcard_constraints:
        count_correction_id=r"ovrlpy"

use rule adaptCorrectedCountsBySplit2Seurat as adaptCorrectedCountsByResolviUnsupervised2Seurat with:
    output:
        protected(f'{config["output_path"]}/post_count_correction_std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{count_correction_id}}/mixture_k={config["count_correction"]["resolvi"]["train"]["mixture_k"]}/num_samples={config["count_correction"]["resolvi"]["predict"]["num_samples"]}/raw_seurat.rds')
    params:
        spatial_dimname=sec.SEURAT_SPATIAL_DIM_NAME,
        segmentation_id=lambda wildcards: wildcards.segmentation_id,
        segmentation_method=lambda wildcards: extract_layers_from_experiments(
            wildcards.segmentation_id,
            0,
            sep_in="_",
        )[0],
        sample_id=lambda wildcards: wildcards.sample_id,
        condition=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            0,
        )[0],
        gene_panel=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            1,
        )[0],
        donor=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            2,
        )[0],
        sample=lambda wildcards: extract_layers_from_experiments(
            wildcards.sample_id,
            3,
        )[0],
        normalisation_id=None,
        annotation_id=None,
        annotation_approach=None,
        annotation_reference_name=None,
        annotation_method=None,
        annotation_level=None,
        annotation_mode=None,
        count_correction_id=lambda wildcards: wildcards.count_correction_id
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
