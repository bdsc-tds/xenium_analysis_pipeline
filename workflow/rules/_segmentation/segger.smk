#######################################
#              Functions              #
#######################################

def get_model_version4runSeggerPredict(wildcards) -> int:
    model_dir: str = f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/trained_model'

    for dir_path, dir_names, _ in os.walk(model_dir):
        if not all(re.match(r"version_\d+", i) for i in dir_names):
            continue
        
        return max(
            [
                int(re.match(r"version_(\d+)", i).group(1))
                for i in dir_names
            ]
        )
    
    return -1

def get_input2_or_params4normaliseSegger(wildcards) -> dict[str, str]:
    ret: dict[str, str] = {
        "data_dir": get_input2_or_params4run10x(wildcards),
        "segmentation": f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/processed_results/segmentation.csv'
    }

    meets_min: bool = get_xeniumranger_version(
        checkpoints.check10xVersions.get(sample_id=wildcards.sample_id).output[0],
        min_version=(10,)
    )

    if meets_min:
        ret["polygons"] = f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/processed_results/segmentation_polygons_feat_col.json'
    else:
        ret["polygons"] = f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/processed_results/segmentation_polygons_geom_col.json'

    return ret


#######################################
#                Rules                #
#######################################

rule runSeggerPreprocess:
    input:
        get_input2_or_params4run10x
    output:
        directory(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data/data')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/runSeggerPreprocess.log'
    params:
        input=lambda wildcards: (
            get_input2_or_params4run10x(
                wildcards,
                for_input=False
            )
        ),
        tile_width=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "preprocess",
            "tile_width",
            replace_none=200
        ),
        tile_height=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "preprocess",
            "tile_height",
            replace_none=200
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "preprocess",
            "_other_options",
            replace_none=""
        )
    threads:
        lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "preprocess",
            "_threads",
            replace_none=1
        )
    resources:
        mem_mb=lambda wildcards, threads, attempt: threads * 4096 * attempt
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n segger_cuda python3 /opt/segger_dev/src/segger/cli/create_dataset_fast.py "
        "--base_dir {params.input} "
        "--data_dir {output} "
        "--sample_type xenium "
        "--tile_width {params.tile_width} "
        "--tile_height {params.tile_height} "
        "--n_workers {threads} "
        "{params.other_options} &> {log}"

rule runSeggerTrain:
    input:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data/data'
    output:
        directory(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/trained_model')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/runSeggerTrain.log'
    params:
        num_tx_tokens=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "num_tx_tokens",
            replace_none=500
        ),
        accelerator=lambda wildcards: "cuda" if _use_gpu() else "cpu",
        devices=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "devices",
            replace_none=0
        ),
        max_epochs=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "max_epochs",
            replace_none=200
        ),
        batch_size=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "batch_size",
            replace_none=4
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "_other_options",
            replace_none=""
        )
    threads:
        lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "train",
            "_threads",
            replace_none=2
        )
    resources:
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        mem_mb=lambda wildcards, threads, attempt: threads * attempt * (2048 if _use_gpu() else 20480),
        slurm_extra=lambda wildcards: get_slurm_extra(
            wildcards,
            get_dict_value(
                config,
                "segmentation",
                "segger",
                "train",
                "devices",
                replace_none=1,
            ),
        )
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n segger_cuda python3 /opt/segger_dev/src/segger/cli/train_model.py "
        "--dataset_dir {input} "
        "--models_dir {output} "
        "--sample_tag run "
        "--num_tx_tokens {params.num_tx_tokens} "
        "--accelerator {params.accelerator} "
        "--devices {params.devices} "
        "--max_epochs {params.max_epochs} "
        "--batch_size {params.batch_size} "
        "--num_workers {threads} "
        "{params.other_options} &> {log}"

rule runSeggerPredict:
    input:
        processed_data=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data/data',
        trained_model=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/trained_model',
        transcripts_file=get_input2_or_params4runProseg
    output:
        directory(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/runSeggerPredict.log'
    params:
        model_version=get_model_version4runSeggerPredict,
        transcripts_file=lambda wildcards: get_input2_or_params4runProseg(
            wildcards,
            for_input=False
        ),
        use_cc=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "predict",
            "use_cc",
            replace_none=False
        ),
        other_options=lambda wildcards: get_dict_value(
            config,
            "segmentation",
            "segger",
            "predict",
            "_other_options",
            replace_none=""
        )
    threads:
        1
    resources:
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        mem_mb=lambda wildcards, attempt: min(
            get_size(
                get_input2_or_params4runProseg(
                    wildcards,
                    for_input=False
                    )
            ) * 10**-6 * attempt * (100 if _use_gpu() else 500),
            1024000
        ),
        slurm_extra=get_slurm_extra
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n segger_cuda python3 /opt/segger_dev/src/segger/cli/predict_fast.py "
        "--segger_data_dir {input.processed_data} "
        "--models_dir {input.trained_model} "
        "--benchmarks_dir {output} "
        "--transcripts_file {params.transcripts_file} "
        "--batch_size 1 "
        "--num_workers {threads} "
        "--model_version {params.model_version} "
        "--knn_method kd_tree "
        "--cell_id_col segger_cell_id "
        "--use_cc {params.use_cc} "
        "{params.other_options} &> {log}"

rule cleanSeggerPredictDir:
    input:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results'
    output:
        protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results/segger_adata.h5ad'),
        protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results/segger_transcripts.parquet'),
        protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results/segmentation_log.json'),
        protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results/compressed_transcripts_df.parquet')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/cleanSeggerPredictDir.log'
    conda:
        "../../envs/pyarrow.yml"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_size(
            input[0]
        ) * 4 * 10**-6 * attempt
    shell:
        "python3 workflow/scripts/_segmentation/clean_segger_predict_results.py "
        "--dir {input} "
        "-l {log}"

rule runSegger2Baysor:
    input:
        data_file=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results/segger_transcripts.parquet',
        xr_version=f'{config["output_path"]}/reprocessed/{{sample_id}}/versions.json'
    output:
        segmentation=protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/processed_results/segmentation.csv'),
        polygons_feat=protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/processed_results/segmentation_polygons_feat_col.json'),
        polygons_geom=protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/processed_results/segmentation_polygons_geom_col.json')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/runSegger2Baysor.log'
    params:
        other_options=lambda wildcards, input: "--no-prior2baysor07" if get_xeniumranger_version(
                input.xr_version,
                min_version=(3, 1)
            ) else "--prior2baysor07"
    resources:
        mem_mb=lambda wildcards, input, attempt: input.size_mb * 60 * attempt
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n segger_cuda python3 workflow/scripts/_segmentation/convert_segger2baysor.py "
        "--inseg {input.data_file} "
        "--outseg {output.segmentation} "
        "--outpolyfeat {output.polygons_feat} "
        "--outpolygeom {output.polygons_geom} "
        "-l {log} "
        "{params.other_options}"

rule normaliseSegger:
    input:
        unpack(get_input2_or_params4normaliseSegger)
    output:
        directory(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/normalised_results')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/normaliseSegger.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/segger/{{sample_id}}',
        abs_input_data_dir=lambda wildcards: os.path.abspath(
            get_input2_or_params4run10x(wildcards, for_input=False)
        ),
        abs_input_segmentation=lambda wildcards: os.path.abspath(
            get_input2_or_params4normaliseSegger(wildcards)["segmentation"]
        ),
        abs_input_polygons=lambda wildcards: os.path.abspath(
            get_input2_or_params4normaliseSegger(wildcards)["polygons"]
        ),
        abs_log=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/logs/normaliseSegger.log'
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

rule zipSeggerPreprocessed:
    input:
        tiles=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data/data',
        pred=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/raw_results'
    output:
        protected(f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data/tiles.tgz')
    log:
        f'{config["output_path"]}/segmentation/segger/{{sample_id}}/logs/zipSeggerPreprocess.log'
    params:
        work_dir=f'{config["output_path"]}/segmentation/segger/{{sample_id}}/preprocessed_data',
        abs_output=lambda wildcards: os.path.abspath(
            f'{config["output_path"]}/segmentation/segger/{wildcards.sample_id}/preprocessed_data/tiles.tgz'
        )
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * 2 * attempt, 2048)
    shell:
        "cd {params.work_dir} && "
        "[[ -f tiles.tgz ]] && rm -f tiles.tgz; "
        "tar --remove-files -czf {params.abs_output} data &> {log}"
