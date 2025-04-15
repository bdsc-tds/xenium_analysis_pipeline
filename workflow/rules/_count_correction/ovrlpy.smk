#######################################
#              Functions              #
#######################################

def get_input2_or_params4runOvrlpy(wildcards, for_input: bool = True) -> str:
    if re.match(
        r"^_proseg$",
        wildcards.segmentation_id4ovrlpy,
        flags=re.IGNORECASE,
    ) is not None:
        return f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/raw_results/transcript-metadata.csv.gz'
    
    if re.match(
        r"^_general$",
        wildcards.segmentation_id4ovrlpy,
        flags=re.IGNORECASE,
    ) is not None:
        use_raw_data, ret = get_raw_data_dir(wildcards.sample_id)

        if for_input:
            return ret

        if use_raw_data:
            return normalise_path(
                ret,
                pat_anchor_file=r"transcripts.parquet",
                pat_flags=re.IGNORECASE,
                return_dir=False,
                check_exist=True,
            )

        return normalise_path(
            ret,
            candidate_paths=("outs",),
            pat_anchor_file=r"transcripts.parquet",
            pat_flags=re.IGNORECASE,
            return_dir=False,
            check_exist=False,
        )

    raise RuntimeError(f'Error! Unknown segmentation id: {wildcards.segmentation_id4ovrlpy}.')


def get_input2_or_params4getCorrectedCountsFromOvrlpy(wildcards, for_input: bool = True) -> dict[str, str]:
    prefix_ovrlpy: str = f'{config["output_path"]}/count_correction'

    ret: dict[str, str] = {}

    segmentation_id4ovrlpy: str = ""

    if re.match(
        r"^proseg_expected$",
        wildcards.segmentation_id,
        flags=re.IGNORECASE,
    ) is not None:
        segmentation_id4ovrlpy = "_proseg"

        ret["transcripts"] = f'{config["output_path"]}/segmentation/proseg/{wildcards.sample_id}/raw_results/transcript-metadata.csv.gz'
    else:
        segmentation_id4ovrlpy = "_general"

        use_raw_data, _ret = get_raw_data_dir(wildcards.sample_id)

        if for_input:
            ret["transcripts"] = _ret
        else:
            if use_raw_data:
                ret["transcripts"] = normalise_path(
                    _ret,
                    pat_anchor_file=r"transcripts.parquet",
                    pat_flags=re.IGNORECASE,
                    return_dir=False,
                    check_exist=True,
                )
            else:
                ret["transcripts"] = normalise_path(
                    _ret,
                    candidate_paths=("outs",),
                    pat_anchor_file=r"transcripts.parquet",
                    pat_flags=re.IGNORECASE,
                    return_dir=False,
                    check_exist=False,
                )

    ret["signal_integrity"] = os.path.join(
        prefix_ovrlpy,
        segmentation_id4ovrlpy,
        wildcards.sample_id,
        "ovrlpy/signal_integrity.parquet",
    )

    ret["transcript_info"] = os.path.join(
        prefix_ovrlpy,
        segmentation_id4ovrlpy,
        wildcards.sample_id,
        "ovrlpy/transcript_info.parquet",
    )

    return ret


#######################################
#                Rules                #
#######################################

rule runOvrlpy:
    input:
        get_input2_or_params4runOvrlpy
    output:
        signal_integrity=protected(f'{config["output_path"]}/count_correction/{{segmentation_id4ovrlpy}}/{{sample_id}}/ovrlpy/signal_integrity.parquet'),
        signal_strength=protected(f'{config["output_path"]}/count_correction/{{segmentation_id4ovrlpy}}/{{sample_id}}/ovrlpy/signal_strength.parquet'),
        transcript_info=protected(f'{config["output_path"]}/count_correction/{{segmentation_id4ovrlpy}}/{{sample_id}}/ovrlpy/transcript_info.parquet')
    params:
        input_transcripts=lambda wildcards: get_input2_or_params4runOvrlpy(
            wildcards,
            for_input=False,
        ),
        proseg_format=lambda wildcards: '--proseg_format' if wildcards.segmentation_id4ovrlpy == '_proseg' else ''
    log:
        f'{config["output_path"]}/count_correction/{{segmentation_id4ovrlpy}}/{{sample_id}}/ovrlpy/logs/runOvrlpy.log'
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
        unpack(get_input2_or_params4getCorrectedCountsFromOvrlpy)
    output:
        corrected_counts=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/corrected_counts.h5'),
        cells_mean_integrity_unfiltered=protected(f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/cells_mean_integrity_unfiltered.parquet')
    params:
        input_transcripts=lambda wildcards: get_dict_value(
            get_input2_or_params4getCorrectedCountsFromOvrlpy(
                wildcards,
                for_input=False,
            ),
            "transcripts",
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
        f'{config["output_path"]}/count_correction/{{segmentation_id}}/{{sample_id}}/ovrlpy/signal_integrity_threshold={config["count_correction"]["ovrlpy"]["signal_integrity_threshold"]}/logs/getCorrectedCountsFromOvrlpy.log'
    container:
        config["containers"]["python_cuda"]
    resources:
        mem_mb=lambda wildcards, attempt: min(
            sum(
                [
                    get_size(i)
                    for i in get_input2_or_params4getCorrectedCountsFromOvrlpy(
                        wildcards,
                        for_input=False,
                    ).values()
                ]
            ) * 1e-6 * attempt,
            512000
        )
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_count_correction/ovrlpy_sample_correction.py "
        "--sample_transcripts_path {params.input_transcripts} "
        "--sample_signal_integrity {input.signal_integrity} "
        "--sample_transcript_info {input.transcript_info} "
        "--out_file_corrected_counts {output.corrected_counts} "
        "--out_file_cells_mean_integrity_unfiltered {output.cells_mean_integrity_unfiltered} "
        "--signal_integrity_threshold {params.signal_integrity_threshold} "
        "{params.proseg_format} "
        "-l {log}"
