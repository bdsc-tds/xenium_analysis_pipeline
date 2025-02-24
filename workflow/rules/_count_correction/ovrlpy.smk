#######################################
#              Functions              #
#######################################

def get_input2_or_params4runOvrlpy(wildcards, for_input: bool = True) -> str:
    prefix: str = f'{config["output_path"]}/segmentation'
    ret: str = ""

    if re.match(
        r"^10x_\w*?_?0um$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            f"10x_0um/{wildcards.sample_id}/normalised_results",
        )

        if not for_input:
            ret = normalise_path(
                ret,
                candidate_paths=("outs",),
                pat_anchor_file=r"transcripts.parquet",
                pat_flags=re.IGNORECASE,
                return_dir=False,
                check_exist=False,
            )
    elif re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        ret = os.path.join(
            prefix,
            f'proseg/{wildcards.sample_id}/raw_results/transcript-metadata.csv.gz',
        )
    else:
        raise RuntimeError(f'Error! Do not run Ovrlpy on samples segmented by: {wildcards.segmentation_id}, except for "10x_0um" and "proseg_expected".')

    return ret


#######################################
#                Rules                #
#######################################

rule runOvrlpy:
    input:
        get_input2_or_params4runOvrlpy
    output:
        signal_integrity=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity.parquet'),
        signal_strength=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_strength.parquet'),
        transcript_info=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/transcript_info.parquet')
    params:
        input_transcripts=lambda wildcards: get_input2_or_params4runOvrlpy(
            wildcards,
            for_input=False,
        ),
        proseg_format=lambda wildcards: '--proseg_format' if wildcards.segmentation_id == 'proseg_expected' else ''
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/logs/runOvrlpy.log'
    wildcard_constraints:
        segmentation_id=r"(10x_\w*?_?0um)|(proseg_expected)"
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_size(
                get_input2_or_params4runOvrlpy(
                    wildcards,
                    for_input=False,
                )
            ) * 1e-6 * attempt**2 * 80,
            1024000
        )
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/ovrlpy_sample.py "
        "--sample_transcripts_path {params.input_transcripts} "
        "--out_file_signal_integrity {output.signal_integrity} "
        "--out_file_signal_strength {output.signal_strength} "
        "--out_file_transcript_info {output.transcript_info} "
        "{params.proseg_format} "
        "-l {log}"

rule getCorrectedCountsFromOvrlpy:
    input:
        transcripts=get_input2_or_params4runOvrlpy,
        signal_integrity=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity.parquet',
        transcript_info=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/transcript_info.parquet'
    output:
        corrected_counts=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/corrected_counts.h5'),
        cells_mean_integrity_filtered=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/cells_mean_integrity_filtered.parquet')
    params:
        input_transcripts=lambda wildcards: get_input2_or_params4runOvrlpy(
            wildcards,
            for_input=False,
        ),
        signal_integrity_threshold=get_dict_value(
            config,
            "count_correction",
            "ovrlpy",
            "signal_integrity_threshold",
            replace_none=0.5,
        ),
        proseg_format=lambda wildcards: '--proseg_format' if wildcards.segmentation_id == 'proseg_expected' else ''
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/logs/getCorrectedCountsFromOvrlpy.log'
    wildcard_constraints:
        segmentation_id=r"(10x_\w*?_?0um)|(proseg_expected)"
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_size(
                get_input2_or_params4runOvrlpy(
                    wildcards,
                    for_input=False,
                )
            ) * 1e-6 * attempt,
            512000
        )
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/ovrlpy_sample_correction.py "
        "--sample_transcripts_path {params.input_transcripts} "
        "--sample_signal_integrity {input.signal_integrity} "
        "--sample_transcript_info {input.transcript_info} "
        "--out_file_corrected_counts {output.corrected_counts} "
        "--out_file_cells_mean_integrity_filtered {output.cells_mean_integrity_filtered} "
        "--signal_integrity_threshold {params.signal_integrity_threshold} "
        "{params.proseg_format} "
        "-l {log}"

rule getUnfilteredCellMeanIntegrityFromOvrlpy:
    input:
        transcripts=get_input2_or_params4runOvrlpy,
        signal_integrity=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity.parquet',
        transcript_info=f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/transcript_info.parquet'
    output:
        protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/cells_mean_integrity_unfiltered.parquet')
    params:
        input_transcripts=lambda wildcards: get_input2_or_params4runOvrlpy(
            wildcards,
            for_input=False,
        ),
        proseg_format=lambda wildcards: '--proseg_format' if wildcards.segmentation_id == 'proseg' else ''
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/logs/getCorrectedCountsFromOvrlpy.log'
    wildcard_constraints:
        segmentation_id=r"(10x_\w*?_?0um)|(proseg_expected)"
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_size(
                get_input2_or_params4runOvrlpy(
                    wildcards,
                    for_input=False,
                )
            ) * 1e-6 * attempt,
            512000
        )
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/ovrlpy_sample_cell_mean_integrity.py "
        "--sample_transcripts_path {params.input_transcripts} "
        "--sample_signal_integrity {input.signal_integrity} "
        "--sample_transcript_info {input.transcript_info} "
        "--out_file_cells_mean_integrity_unfiltered {output} "
        "{params.proseg_format} "
        "-l {log}"