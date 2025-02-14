from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# params
num_samples = 100

out_files = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    if segmentation.stem == 'proseg':
                        path = sample / 'raw_results'
                        segmentation_name = 'proseg_expected'
                    else:
                        path = sample / "normalised_results/outs"
                        segmentation_name = segmentation.stem
                    
                    k = (segmentation_name,condition.stem,panel.stem,donor.stem,sample.stem)
                    name = '/'.join(k)


                    if path.exists():

                        out_file_resolvi_corrected_counts = results_dir / f'resolvi/{name}/resolvi_corrected_counts.h5'
                        out_file_resolvi_proportions = results_dir/ f'resolvi/{name}/resolvi_proportions.parquet'
                        out_files.extend([out_file_resolvi_corrected_counts,out_file_resolvi_proportions])

                        rule:
                            name: f'resolvi/{name}'
                            input:
                                path=path,
                            output:
                                out_file_resolvi_corrected_counts=out_file_resolvi_corrected_counts,
                                out_file_resolvi_proportions=out_file_resolvi_proportions,
                            params:
                                min_counts=min_counts,
                                min_features=min_features,
                                max_counts=max_counts,
                                max_features=max_features,
                                min_cells=min_cells,
                                num_samples=num_samples,
                            threads: 1
                            resources:
                                mem='400GB',
                                runtime='8h',
                                slurm_partition = "gpu",
                                slurm_extra = '--gres=gpu:1',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file_resolvi_corrected_counts})"

                                python workflow/scripts/xenium/resolvi_sample.py \
                                --path {input.path} \
                                --out_file_resolvi_corrected_counts {output.out_file_resolvi_corrected_counts} \
                                --out_file_resolvi_proportions {output.out_file_resolvi_proportions} \
                                --min_counts {params.min_counts} \
                                --min_features {params.min_features} \
                                --max_counts {params.max_counts} \
                                --max_features {params.max_features} \
                                --min_cells {params.min_cells} \
                                --num_samples {params.num_samples} \

                                echo "DONE"
                                """


rule resolvi_all:
    input:
        out_files