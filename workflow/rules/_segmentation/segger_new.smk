#######################################
#              Overview               #
#######################################
# The *new* Segger (dpeerlab/segger) segmentation arm ("segger_new").
#
# Unlike the old Segger (bdsc-tds fork, 5 rules: preprocess -> train -> predict
# -> clean -> convert, baked into skang's python_cuda.sif), the new Segger is a
# single unified command `segger segment -i <xenium_dir> -o <out>` that does
# preprocess+train+predict in one Lightning run. We run it from a host conda env
# (workflow/envs/segger_new_cuda.yml) via Snakemake's `conda:` directive instead
# of a container.
#
# Chain (mirrors the other methods and terminates in the SAME xeniumranger
# import-segmentation step, so `normalised_results/outs/` is identical to every
# other method):
#   runSeggerNew        : segger segment  -> raw_results/segger_segmentation.parquet
#   runSeggerNew2Baysor : reshape output  -> processed_results/segmentation.csv + polygon JSONs
#   normaliseSeggerNew  : xeniumranger import-segmentation -> normalised_results/
#
# Reused helpers (defined in earlier-included .smk files, same namespace):
#   get_input2_or_params4run10x     (10x.smk)     -> raw Xenium bundle dir (auto-detected by segger)
#   get_input2_or_params4runProseg  (proseg.smk)  -> raw transcripts.parquet path (for row_index join)
#   get_xeniumranger_version        (baysor.smk)
#   _use_gpu / get_slurm_gpu_partition_name / get_slurm_extra / get_dict_value (segmentation.smk)
#   get_size (run_time_utils.py)

#######################################
#              Functions              #
#######################################

def get_input2_or_params4normaliseSeggerNew(wildcards) -> dict[str, str]:
    ret: dict[str, str] = {
        "data_dir": get_input2_or_params4run10x(wildcards),
        "segmentation": f'{config["output_path"]}/segmentation/segger_new/{wildcards.sample_id}/processed_results/segmentation.csv'
    }

    meets_min: bool = get_xeniumranger_version(
        checkpoints.check10xVersions.get(sample_id=wildcards.sample_id).output[0],
        min_version=(10,)
    )

    if meets_min:
        ret["polygons"] = f'{config["output_path"]}/segmentation/segger_new/{wildcards.sample_id}/processed_results/segmentation_polygons_feat_col.json'
    else:
        ret["polygons"] = f'{config["output_path"]}/segmentation/segger_new/{wildcards.sample_id}/processed_results/segmentation_polygons_geom_col.json'

    return ret


#######################################
#                Rules                #
#######################################

rule runSeggerNew:
    input:
        get_input2_or_params4run10x
    output:
        # Declare the parquet itself, NOT directory(raw_results): Snakemake does
        # not resolve a file *inside* another rule's directory() output, so
        # runSeggerNew2Baysor would fail with MissingInputException on a fresh
        # run. Same reason the old segger arm has cleanSeggerPredictResults.
        # (`lightning_logs/` is also written into raw_results but is untracked.)
        # Deliberately not protected(): this arm exists to be re-run across
        # n_epochs / cells_representation / segmentation_loss settings, and a
        # write-protected output makes every such re-run need a manual chmod.
        f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/raw_results/segger_segmentation.parquet'
    log:
        f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/logs/runSeggerNew.log'
    params:
        out_dir=f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/raw_results',
        input_dir=lambda wildcards: get_input2_or_params4run10x(
            wildcards,
            for_input=False
        ),
        n_epochs=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger_new",
            "n_epochs",
            replace_none=20
        ),
        cells_representation=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger_new",
            "cells_representation",
            replace_none="pca"
        ),
        segmentation_loss=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger_new",
            "segmentation_loss",
            replace_none="triplet"
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger_new",
            "_other_options",
            replace_none=""
        )
    threads:
        lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger_new",
            "_threads",
            replace_none=4
        )
    resources:
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        # Scaling copied from runSeggerPredict, plus a floor (cf. runBaysor2Xenium).
        # Peak RSS over a 36-sample cohort tracked the transcript-parquet size at
        # <= 93 MB per MB, so the 100x term holds, but the ratio climbs as samples
        # shrink -- there is a fixed overhead the linear term does not model. The
        # smallest sample (23 MB parquet) is sized at 2.2 GB by the formula alone
        # and is OOM-killed at ~1.9 GB; comparable samples peaked at 4-5 GB.
        mem_mb=lambda wildcards, attempt: max(
            min(
                get_size(
                    get_input2_or_params4runProseg(
                        wildcards,
                        for_input=False
                        )
                ) * 10**-6 * attempt * (100 if _use_gpu() else 500),
                1024000
            ),
            8192 * attempt
        ),
        slurm_extra=get_slurm_extra
    conda:
        "../../envs/segger_new_cuda.yml"
    shell:
        # segger's CLI validates that -o already exists.
        "mkdir -p {params.out_dir} && "
        # Fix CUDA lib discovery for the pip RAPIDS/cupy stack in the bare conda
        # env (LD_LIBRARY_PATH for cupy's libcudart.so.12 + unversioned symlinks
        # for cudf/numba). See the sourced helper for details.
        "source workflow/scripts/_segmentation/segger_cuda_env.sh && "
        "LOG_LEVEL=INFO segger segment "
        "-i {params.input_dir} "
        "-o {params.out_dir} "
        "--n-epochs {params.n_epochs} "
        "--cells-representation {params.cells_representation} "
        "--segmentation-loss {params.segmentation_loss} "
        # segger's own default is save_anndata=True (from ISTSegmentationWriter);
        # the downstream converter only reads segger_segmentation.parquet, so
        # skip the extra cell x gene h5ad write. This matches how the 36-sample
        # benchmark run was executed.
        "--no-save-anndata "
        "{params.other_options} &> {log}"

rule runSeggerNew2Baysor:
    input:
        data_file=f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/raw_results/segger_segmentation.parquet',
        xr_version=f'{config["output_path"]}/reprocessed/{{sample_id}}/versions.json'
    output:
        segmentation=protected(f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/processed_results/segmentation.csv'),
        polygons_feat=protected(f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/processed_results/segmentation_polygons_feat_col.json'),
        polygons_geom=protected(f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/processed_results/segmentation_polygons_geom_col.json')
    log:
        f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/logs/runSeggerNew2Baysor.log'
    params:
        # raw transcripts.parquet -> recover transcript_id / cell_id / x / y via row_index
        transcripts_file=lambda wildcards: get_input2_or_params4runProseg(
            wildcards,
            for_input=False
        ),
        other_options=lambda wildcards, input: "--no-prior2baysor07" if get_xeniumranger_version(
                input.xr_version,
                min_version=(3, 1)
            ) else "--prior2baysor07"
    resources:
        # NOT sized off input.size_mb (the runSegger2Baysor idiom): that rule's
        # input is the old segger_transcripts.parquet, which carries the
        # transcript ids and coordinates and is ~5.6x larger than the new, lean
        # segger_segmentation.parquet (586 MB vs 105 MB on the same sample). The
        # cost here is dominated by the join against the raw transcripts and the
        # per-cell boundary generation, so size off the raw transcripts (a param,
        # hence invisible to input.size_mb) and keep a floor.
        mem_mb=lambda wildcards, attempt: max(
            get_size(
                get_input2_or_params4runProseg(
                    wildcards,
                    for_input=False
                    )
            ) * 10**-6 * attempt * 30,
            8192 * attempt
        )
    conda:
        "../../envs/segger_new_cuda.yml"
    shell:
        "python3 workflow/scripts/_segmentation/convert_seggernew2baysor.py "
        "--inseg {input.data_file} "
        "--intx {params.transcripts_file} "
        "--outseg {output.segmentation} "
        "--outpolyfeat {output.polygons_feat} "
        "--outpolygeom {output.polygons_geom} "
        "-l {log} "
        "{params.other_options}"

rule normaliseSeggerNew:
    input:
        unpack(get_input2_or_params4normaliseSeggerNew)
    output:
        directory(f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}/logs/normaliseSeggerNew.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/segger_new/{{sample_id}}',
        abs_input_data_dir=lambda wildcards: os.path.abspath(
            get_input2_or_params4run10x(wildcards, for_input=False)
        ),
        abs_input_segmentation=lambda wildcards: os.path.abspath(
            get_input2_or_params4normaliseSeggerNew(wildcards)["segmentation"]
        ),
        abs_input_polygons=lambda wildcards: os.path.abspath(
            get_input2_or_params4normaliseSeggerNew(wildcards)["polygons"]
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/segger_new/{wildcards.sample_id}/logs/normaliseSeggerNew.log'
        ),
        localmem=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_memory"
        )
    retries:
        0
    threads:
        get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_threads"
        )
    resources:
        mem_mb=get_dict_value(
            config,
            "segmentation",
            "_normalisation",
            "_memory"
        ) * 1024
    container:
        config["containers"]["10x"]
    shell:
        "cd {params.work_dir} && "
        "xeniumranger import-segmentation --id=normalised_results "
        "--xenium-bundle {params.abs_input_data_dir} "
        "--transcript-assignment={params.abs_input_segmentation} "
        "--viz-polygons={params.abs_input_polygons} "
        "--units=microns "
        "--localcores={threads} "
        "--localmem={params.localmem} &> {params.abs_log}"
