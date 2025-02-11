from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
signal_integrity_threshold = 0.5
ref_segmentation = xenium_dir / '10x_0um' # ovrlpy output does not depend on segmentation, just run for 10x_0um

out_files_ovrlpy = []
for condition in (conditions := ref_segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):

                k = (ref_segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                sample_transcripts_path = sample / "normalised_results/outs/transcripts.parquet"
                name = '/'.join(k)

                if sample_transcripts_path.exists():

                    out_file_signal_integrity = results_dir / f'ovrlpy/{name}/signal_integrity.parquet' 
                    out_file_signal_strength = results_dir / f'ovrlpy/{name}/signal_strength.parquet'
                    out_file_transcript_info = results_dir / f'ovrlpy/{name}/transcript_info.parquet'
                    # out_file_doublet_df = results_dir / f'ovrlpy/{name}/doublet_df.parquet'

                    out_files_ovrlpy.extend([
                                    out_file_signal_integrity,
                                    out_file_signal_strength,
                                    out_file_transcript_info,
                                        #out_file_doublet_df
                                        ])

                    rule:
                        name: f'ovrlpy/{name}'
                        input:
                            sample_transcripts_path=sample_transcripts_path,
                        output:
                            out_file_signal_integrity=out_file_signal_integrity,
                            out_file_signal_strength=out_file_signal_strength,
                            out_file_transcript_info=out_file_transcript_info,
                            # out_file_doublet_df=out_file_doublet_df,
                        threads: 1
                        resources:
                            mem='400GB',
                            runtime='15h',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_signal_integrity})"

                            python workflow/scripts/xenium/ovrlpy_sample.py \
                            --sample_transcripts_path {input.sample_transcripts_path} \
                            --out_file_signal_integrity {output.out_file_signal_integrity} \
                            --out_file_signal_strength {output.out_file_signal_strength} \
                            --out_file_transcript_info {output.out_file_transcript_info} \

                            echo "DONE"
                            """
                            # --out_file_doublet_df {output.out_file_doublet_df} \


out_files_ovrlpy_correction = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    if sample.stem == '1FYB':
                        continue

                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                    k_ref = (ref_segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)

                    name = '/'.join(k)
                    ref_name = '/'.join(k_ref)

                    sample_transcripts_path = sample / "normalised_results/outs/transcripts.parquet"
                    sample_signal_integrity = results_dir / f'ovrlpy/{ref_name}/signal_integrity.parquet' 
                    sample_transcript_info = results_dir / f'ovrlpy/{ref_name}/transcript_info.parquet'

                    if sample_transcripts_path.exists():
                        out_file_corrected_counts = results_dir / f'ovrlpy_correction/{name}/corrected_counts_{signal_integrity_threshold=}.h5'
                        out_file_cells_mean_integrity = results_dir / f'ovrlpy_correction/{name}/cells_mean_integrity.parquet'

                        out_files_ovrlpy_correction.extend([ out_file_corrected_counts, out_file_cells_mean_integrity,])

                        rule:
                            name: f'ovrlpy_correction/{name}'
                            input:
                                sample_transcripts_path=sample_transcripts_path,
                                sample_signal_integrity=sample_signal_integrity,
                                sample_transcript_info=sample_transcript_info,
                            output:
                                out_file_corrected_counts=out_file_corrected_counts,
                                out_file_cells_mean_integrity=out_file_cells_mean_integrity,
                            params:
                                signal_integrity_threshold=signal_integrity_threshold,
                            threads: 1
                            resources:
                                mem='30GB' if panel.stem == '5k' else '20GB',
                                runtime='20m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file_corrected_counts})"

                                python workflow/scripts/xenium/ovrlpy_sample_correction.py \
                                --sample_transcripts_path {input.sample_transcripts_path} \
                                --sample_signal_integrity {input.sample_signal_integrity} \
                                --sample_transcript_info {input.sample_transcript_info} \
                                --out_file_corrected_counts {output.out_file_corrected_counts} \
                                --out_file_cells_mean_integrity {output.out_file_cells_mean_integrity} \
                                --signal_integrity_threshold {params.signal_integrity_threshold} \

                                echo "DONE"
                                """



rule ovrlpy_all:
    input:
        out_files_ovrlpy

rule ovrlpy_correction_all:
    input:
        out_files_ovrlpy_correction