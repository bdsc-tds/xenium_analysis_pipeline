# Python implementation of RCTD (rctd-py), included only when
# `cell_type_annotation.reference_based.rctd.language` is `python`.
#
# The chain writes the same files, at the same paths, as `runReferenceBasedRCTD` (rctd.smk), so
# post-processing and every downstream rule run unchanged:
#   1. runConvertSeuratReferenceToH5: Seurat .rds reference -> {stem}_counts.h5 + {stem}_metadata.parquet,
#      next to the .rds; only scheduled when the configured reference path ends in .rds. With sctransform,
#      runConvertSeuratReferenceSCTToH5 also writes the SCT assay counts to {stem}_SCT_counts.h5.
#   2. exportQuery4runReferenceBasedRCTDPy: Xenium Seurat object -> 10x h5 counts + coordinates.
#   3. runReferenceBasedRCTDPy: rctd-py in doublet mode (GPU, rctd_py env of the python_cuda container),
#      same reference filtering and settings as rctd.R.
#   4. assembleReferenceBasedRCTDPy: spacexr RCTD object + output parquets, as written by rctd.R.
# The RCTD runs after count correction (SPLIT, ovrlpy, resolVI) switch the same way: steps 2-4 run on the
# corrected query, with the query export and the rctd-py outputs in the output folder of the R rule they replace.

#######################################
#              Functions              #
#######################################

def _get_reference_paths4rctd_py() -> list[str]:
    _paths: dict = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_REFERENCE_FILES_NAME,
    )

    return sorted({p for v in _paths.values() for p in v.values()})


def get_reference_files4runReferenceBasedRCTDPy(wildcards) -> dict[str, str]:
    """
    Reference counts and metadata: converted from the .rds when the configured path is a .rds,
    used directly when it is a `{stem}_counts.h5` (labels in `{stem}_metadata.parquet`).
    With sctransform, the counts are the reference's SCT assay, in `{stem}_SCT_counts.h5`.
    """
    path: str = get_path2reference4reference_based_annotation(wildcards)

    if path.endswith(".rds"):
        stem: str = path[: -len(".rds")]
    elif path.endswith("_counts.h5"):
        stem: str = path[: -len("_counts.h5")]
    else:
        raise RuntimeError(f"Error! The reference for rctd-py must be a .rds or a {{stem}}_counts.h5 file: {path}.")

    # Same reference assay as rctd.R: RNA for lognorm, SCT for sctransform
    assay_suffix: str = "" if get_assay_name4annotation(wildcards, True) == cac.REF_SEURAT_DEFAULT_ASSAY else f"_{cac.REF_SEURAT_ALT_ASSAY}"

    return {
        "ref_counts": f"{stem}{assay_suffix}_counts.h5",
        "ref_metadata": f"{stem}_metadata.parquet",
    }


def get_param4runReferenceBasedRCTDPy(wildcards, key: str, default=None):
    return get_dict_value(
        config,
        "cell_type_annotation",
        extract_layers_from_experiments(wildcards.annotation_id, [0])[0],
        "rctd",
        extract_layers_from_experiments(wildcards.annotation_id, [4])[0],
        key,
        replace_none=default,
        inexist_key_ok=True,
    )


def get_donor4runReferenceBasedRCTDPy(wildcards) -> str:
    # sample_id is condition/gene_panel/donor/sample
    return extract_layers_from_experiments(wildcards.sample_id, [2])[0]


def is_int_counts4runReferenceBasedRCTDPy(wildcards) -> bool:
    # Same rule as is_int_counts4runReferenceBasedRCTD for the base (non count-corrected) query
    return re.match("proseg_expected", wildcards.segmentation_id, flags=re.IGNORECASE) is None


# Only the configured .rds references can be converted
_RCTD_PY_RDS_STEMS: list[str] = [p[: -len(".rds")] for p in _get_reference_paths4rctd_py() if p.endswith(".rds")]

if len(_RCTD_PY_RDS_STEMS) > 0:
    wildcard_constraints:
        reference_stem="|".join(re.escape(s) for s in _RCTD_PY_RDS_STEMS)

    include: "../../_reference_preprocessing/convert_seurat_2_h5.smk"


#######################################
#                Rules                #
#######################################

ruleorder: assembleReferenceBasedRCTDPy > runReferenceBasedRCTD

rule exportQuery4runReferenceBasedRCTDPy:
    input:
        query=f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/preprocessed/preprocessed_seurat.rds'
    output:
        counts=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/counts.h5',
        metadata=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/metadata.parquet'
    params:
        xe_assay=lambda wildcards: get_assay_name4annotation(wildcards, False)
    resources:
        mem_mb=lambda wildcards, input, attempt: min(input.size_mb * attempt**2 * 20, 300000)
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/logs/exportQuery4runReferenceBasedRCTDPy.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/rctd_py_export_query.R"


rule runReferenceBasedRCTDPy:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        query_counts=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/counts.h5',
        query_metadata=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/metadata.parquet'
    output:
        results_df=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/results_df.parquet',
        weights=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/weights.parquet',
        weights_doublet=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/weights_doublet.parquet',
        singlet_scores=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/singlet_scores.parquet',
        reference_cells=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/reference_cells.parquet',
        class_df=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/class_df.parquet',
        run_info=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/rctd_py/run_info.json'
    params:
        annotation_level=lambda wildcards: extract_layers_from_experiments(wildcards.annotation_id, [3])[0],
        class_level=lambda wildcards: get_class_level4runReferenceBasedRCTD(wildcards) or "",
        donor=get_donor4runReferenceBasedRCTDPy,
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
        ),
        UMI_min_sigma=lambda wildcards: get_param4runReferenceBasedRCTDPy(wildcards, "UMI_min_sigma"),
        batch_size=lambda wildcards: get_param4runReferenceBasedRCTDPy(wildcards, "_py_batch_size", 10000),
        dtype=lambda wildcards: get_param4runReferenceBasedRCTDPy(wildcards, "_py_dtype", "float64"),
        int_counts_flag=lambda wildcards: "--is_int_counts" if is_int_counts4runReferenceBasedRCTDPy(wildcards) else ""
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd_.+"
    threads:
        lambda wildcards: get_param4runReferenceBasedRCTDPy(wildcards, "_py_threads", 1)
    resources:
        mem_mb=lambda wildcards, input, attempt: min(max(input.size_mb * 10, 16000) * attempt, 500000),
        runtime=lambda wildcards, attempt: 240 * attempt,
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        slurm_extra=get_slurm_extra
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/runReferenceBasedRCTDPy.log'
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run --no-capture-output -n rctd_py python3 workflow/scripts/_cell_type_annotation/_reference_based/rctd_py.py "
        "--ref_counts {input.ref_counts:q} "
        "--ref_metadata {input.ref_metadata:q} "
        "--query_counts {input.query_counts:q} "
        "--query_metadata {input.query_metadata:q} "
        "--annotation_level {params.annotation_level:q} "
        "--class_level={params.class_level:q} "
        "--donor {params.donor:q} "
        "--ref_min_umi {params.REF_MIN_UMI} "
        "--ref_max_umi {params.REF_MAX_UMI} "
        "--xe_min_umi {params.XE_MIN_UMI} "
        "--xe_min_counts {params.XE_MIN_counts} "
        "--cell_min_instance {params.cell_min_instance} "
        "--umi_min_sigma {params.UMI_min_sigma} "
        "--batch_size {params.batch_size} "
        "--dtype {params.dtype} "
        "{params.int_counts_flag} "
        "--out_results_df {output.results_df:q} "
        "--out_weights {output.weights:q} "
        "--out_weights_doublet {output.weights_doublet:q} "
        "--out_singlet_scores {output.singlet_scores:q} "
        "--out_reference_cells {output.reference_cells:q} "
        "--out_class_df {output.class_df:q} "
        "--out_run_info {output.run_info:q} "
        "&> {log:q}"


rule assembleReferenceBasedRCTDPy:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        query_counts=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/counts.h5',
        query_metadata=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/rctd_py_query/metadata.parquet',
        results_df=rules.runReferenceBasedRCTDPy.output.results_df,
        weights=rules.runReferenceBasedRCTDPy.output.weights,
        weights_doublet=rules.runReferenceBasedRCTDPy.output.weights_doublet,
        singlet_scores=rules.runReferenceBasedRCTDPy.output.singlet_scores,
        reference_cells=rules.runReferenceBasedRCTDPy.output.reference_cells,
        class_df=rules.runReferenceBasedRCTDPy.output.class_df,
        run_info=rules.runReferenceBasedRCTDPy.output.run_info
    output:
        rds_output=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output.rds',
        labels=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/labels.parquet',
        scores=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/scores.parquet',
        out_res_df=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/results_df.parquet',
        out_w=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights.parquet',
        out_wd=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/weights_doublet.parquet',
        out_sc=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/singlet_score.parquet',
        out_sm=f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/output/score_mat.parquet'
    params:
        annotation_level=lambda wildcards: extract_layers_from_experiments(wildcards.annotation_id, [3])[0],
        UMI_min_sigma=lambda wildcards: get_param4runReferenceBasedRCTDPy(wildcards, "UMI_min_sigma"),
        cell_min_instance=lambda wildcards: get_dict_value(
            config,
            "experiments",
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
            cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_CELL_MIN_INSTANCES_NAME,
            extract_layers_from_experiments(wildcards.sample_id, [0])[0],
            extract_layers_from_experiments(wildcards.annotation_id, [0, 1])[0],
            replace_none=25,
        ),
        REF_MIN_UMI=cac.REF_MIN_UMI,
        XE_MIN_UMI=cac.XE_MIN_UMI,
        XE_MIN_counts=cac.XE_MIN_counts,
        is_int_counts=is_int_counts4runReferenceBasedRCTDPy
    wildcard_constraints:
        annotation_id=r"reference_based/.+/rctd_.+"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(max(input.size_mb * 20, 32000) * attempt, 300000),
        runtime=lambda wildcards, attempt: 30 * attempt
    log:
        f'{config["output_path"]}/cell_type_annotation/{{segmentation_id}}/{{sample_id}}/{{normalisation_id}}/{{annotation_id}}/logs/assembleReferenceBasedRCTDPy.log'
    container:
        config["containers"]["r"]
    script:
        "../../../scripts/_cell_type_annotation/_reference_based/rctd_py_assemble.R"


#######################################
#   After count correction (Python)   #
#######################################

# File names of the three steps inside an output folder; the counts are not integers after count correction

def _rctd_py_query_files(out_dir: str, prefix: str = "") -> dict[str, str]:
    return {
        f"{prefix}counts": f"{out_dir}/rctd_py_query/counts.h5",
        f"{prefix}metadata": f"{out_dir}/rctd_py_query/metadata.parquet",
    }


def _rctd_py_run_files(out_dir: str) -> dict[str, str]:
    return {
        "results_df": f"{out_dir}/rctd_py/results_df.parquet",
        "weights": f"{out_dir}/rctd_py/weights.parquet",
        "weights_doublet": f"{out_dir}/rctd_py/weights_doublet.parquet",
        "singlet_scores": f"{out_dir}/rctd_py/singlet_scores.parquet",
        "reference_cells": f"{out_dir}/rctd_py/reference_cells.parquet",
        "class_df": f"{out_dir}/rctd_py/class_df.parquet",
        "run_info": f"{out_dir}/rctd_py/run_info.json",
    }


def _rctd_py_output_files(out_dir: str) -> dict[str, str]:
    return {
        "rds_output": f"{out_dir}/output.rds",
        "labels": f"{out_dir}/labels.parquet",
        "scores": f"{out_dir}/scores.parquet",
        "out_res_df": f"{out_dir}/output/results_df.parquet",
        "out_w": f"{out_dir}/output/weights.parquet",
        "out_wd": f"{out_dir}/output/weights_doublet.parquet",
        "out_sc": f"{out_dir}/output/singlet_score.parquet",
        "out_sm": f"{out_dir}/output/score_mat.parquet",
    }

# Split
_RCTD_PY_SPLIT_DIR: str = os.path.dirname(rules.runPostCountCorrectionBySplitReferenceBasedRCTD.output.rds_output)

use rule exportQuery4runReferenceBasedRCTDPy as exportQuery4runPostCountCorrectionBySplitReferenceBasedRCTDPy with:
    input:
        query=lambda wildcards: get_path2query4annotation(wildcards, is_post_correction=True)
    output:
        **_rctd_py_query_files(_RCTD_PY_SPLIT_DIR)
    log:
        f"{_RCTD_PY_SPLIT_DIR}/logs/exportQuery4runPostCountCorrectionBySplitReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTDPy as runPostCountCorrectionBySplitReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_SPLIT_DIR, "query_")
    output:
        **_rctd_py_run_files(_RCTD_PY_SPLIT_DIR)
    params:
        int_counts_flag=""
    log:
        f"{_RCTD_PY_SPLIT_DIR}/logs/runPostCountCorrectionBySplitReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule assembleReferenceBasedRCTDPy as assemblePostCountCorrectionBySplitReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_SPLIT_DIR, "query_"),
        **_rctd_py_run_files(_RCTD_PY_SPLIT_DIR)
    output:
        **_rctd_py_output_files(_RCTD_PY_SPLIT_DIR)
    params:
        is_int_counts=False
    log:
        f"{_RCTD_PY_SPLIT_DIR}/logs/assemblePostCountCorrectionBySplitReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"split.+",
        annotation_id=r"reference_based/.+/rctd_.+"

ruleorder: assemblePostCountCorrectionBySplitReferenceBasedRCTDPy > runPostCountCorrectionBySplitReferenceBasedRCTD

# Ovrlpy
_RCTD_PY_OVRLPY_DIR: str = os.path.dirname(rules.runPostCountCorrectionByOvrlpyReferenceBasedRCTD.output.rds_output)

use rule exportQuery4runReferenceBasedRCTDPy as exportQuery4runPostCountCorrectionByOvrlpyReferenceBasedRCTDPy with:
    input:
        query=lambda wildcards: get_path2query4annotation(wildcards, is_post_correction=True)
    output:
        **_rctd_py_query_files(_RCTD_PY_OVRLPY_DIR)
    log:
        f"{_RCTD_PY_OVRLPY_DIR}/logs/exportQuery4runPostCountCorrectionByOvrlpyReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTDPy as runPostCountCorrectionByOvrlpyReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_OVRLPY_DIR, "query_")
    output:
        **_rctd_py_run_files(_RCTD_PY_OVRLPY_DIR)
    params:
        int_counts_flag=""
    log:
        f"{_RCTD_PY_OVRLPY_DIR}/logs/runPostCountCorrectionByOvrlpyReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule assembleReferenceBasedRCTDPy as assemblePostCountCorrectionByOvrlpyReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_OVRLPY_DIR, "query_"),
        **_rctd_py_run_files(_RCTD_PY_OVRLPY_DIR)
    output:
        **_rctd_py_output_files(_RCTD_PY_OVRLPY_DIR)
    params:
        is_int_counts=False
    log:
        f"{_RCTD_PY_OVRLPY_DIR}/logs/assemblePostCountCorrectionByOvrlpyReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"ovrlpy",
        annotation_id=r"reference_based/.+/rctd_.+"

ruleorder: assemblePostCountCorrectionByOvrlpyReferenceBasedRCTDPy > runPostCountCorrectionByOvrlpyReferenceBasedRCTD

# ResolviUnsupervised
_RCTD_PY_RESOLVIUNSUPERVISED_DIR: str = os.path.dirname(rules.runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTD.output.rds_output)

use rule exportQuery4runReferenceBasedRCTDPy as exportQuery4runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy with:
    input:
        query=lambda wildcards: get_path2query4annotation(wildcards, is_post_correction=True)
    output:
        **_rctd_py_query_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR)
    log:
        f"{_RCTD_PY_RESOLVIUNSUPERVISED_DIR}/logs/exportQuery4runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTDPy as runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR, "query_")
    output:
        **_rctd_py_run_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR)
    params:
        int_counts_flag=""
    log:
        f"{_RCTD_PY_RESOLVIUNSUPERVISED_DIR}/logs/runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule assembleReferenceBasedRCTDPy as assemblePostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR, "query_"),
        **_rctd_py_run_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR)
    output:
        **_rctd_py_output_files(_RCTD_PY_RESOLVIUNSUPERVISED_DIR)
    params:
        is_int_counts=False
    log:
        f"{_RCTD_PY_RESOLVIUNSUPERVISED_DIR}/logs/assemblePostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_unsupervised",
        annotation_id=r"reference_based/.+/rctd_.+"

ruleorder: assemblePostCountCorrectionByResolviUnsupervisedReferenceBasedRCTDPy > runPostCountCorrectionByResolviUnsupervisedReferenceBasedRCTD

# ResolviSupervised
_RCTD_PY_RESOLVISUPERVISED_DIR: str = os.path.dirname(rules.runPostCountCorrectionByResolviSupervisedReferenceBasedRCTD.output.rds_output)

use rule exportQuery4runReferenceBasedRCTDPy as exportQuery4runPostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy with:
    input:
        query=lambda wildcards: get_path2query4annotation(wildcards, is_post_correction=True)
    output:
        **_rctd_py_query_files(_RCTD_PY_RESOLVISUPERVISED_DIR)
    log:
        f"{_RCTD_PY_RESOLVISUPERVISED_DIR}/logs/exportQuery4runPostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule runReferenceBasedRCTDPy as runPostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_RESOLVISUPERVISED_DIR, "query_")
    output:
        **_rctd_py_run_files(_RCTD_PY_RESOLVISUPERVISED_DIR)
    params:
        int_counts_flag=""
    log:
        f"{_RCTD_PY_RESOLVISUPERVISED_DIR}/logs/runPostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r"reference_based/.+/rctd_.+"

use rule assembleReferenceBasedRCTDPy as assemblePostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy with:
    input:
        unpack(get_reference_files4runReferenceBasedRCTDPy),
        **_rctd_py_query_files(_RCTD_PY_RESOLVISUPERVISED_DIR, "query_"),
        **_rctd_py_run_files(_RCTD_PY_RESOLVISUPERVISED_DIR)
    output:
        **_rctd_py_output_files(_RCTD_PY_RESOLVISUPERVISED_DIR)
    params:
        is_int_counts=False
    log:
        f"{_RCTD_PY_RESOLVISUPERVISED_DIR}/logs/assemblePostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy.log"
    wildcard_constraints:
        count_correction_id=r"resolvi_supervised",
        annotation_id=r"reference_based/.+/rctd_.+"

ruleorder: assemblePostCountCorrectionByResolviSupervisedReferenceBasedRCTDPy > runPostCountCorrectionByResolviSupervisedReferenceBasedRCTD
