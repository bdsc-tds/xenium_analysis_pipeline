from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
cell_type_annotation_dir = Path(config['cell_type_annotation_dir'])

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# params
normalisation_method = 'lognorm'
num_samples = 30
batch_size = 1000
macro_batch_size = 50000
mixture_k = 50
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2']

out_files = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    for reference in references:
                        for method in methods:
                            for level in levels:

                                if segmentation.stem == 'proseg':
                                    path = sample / 'raw_results'
                                else:
                                    path = sample / "normalised_results/outs"
                                
                                cell_type_labels=cell_type_annotation_dir / name / f"{normalisation_method}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet"

                                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,reference,method,level)
                                name = '/'.join(k)

                                if path.exists():

                                    out_file_resolvi_corrected_counts = results_dir / f'resolvi_supervised/{name}/corrected_counts.h5'
                                    out_file_resolvi_proportions = results_dir/ f'resolvi_supervised/{name}/proportions.parquet'
                                    out_dir_resolvi_model = results_dir / f'resolvi_supervised/{name}/model/'
                                    out_files.extend([out_file_resolvi_corrected_counts,out_file_resolvi_proportions])

                                    rule:
                                        name: f'resolvi_supervised/{name}'
                                        input:
                                            path=path,
                                        output:
                                            out_file_resolvi_corrected_counts=out_file_resolvi_corrected_counts,
                                            out_file_resolvi_proportions=out_file_resolvi_proportions,
                                            out_dir_resolvi_model=out_dir_resolvi_model,
                                        params:
                                            cell_type_labels=cell_type_labels,
                                            min_counts=min_counts,
                                            min_features=min_features,
                                            max_counts=max_counts,
                                            max_features=max_features,
                                            min_cells=min_cells,
                                            num_samples=num_samples,
                                            batch_size=batch_size,
                                            macro_batch_size=macro_batch_size,
                                            mixture_k=mixture_k,
                                        threads: 1
                                        resources:
                                            mem='80GB',# if panel.stem == '5k' else '10GB',
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
                                            --out_dir_resolvi_model {params.out_dir_resolvi_model}
                                            --cell_type_labels {params.cell_type_labels} \
                                            --min_counts {params.min_counts} \
                                            --min_features {params.min_features} \
                                            --max_counts {params.max_counts} \
                                            --max_features {params.max_features} \
                                            --min_cells {params.min_cells} \
                                            --num_samples {params.num_samples} \
                                            --batch_size {params.batch_size} \
                                            --macro_batch_size {params.macro_batch_size} \
                                            --mixture_k {params.mixture_k}

                                            echo "DONE"
                                            """


rule resolvi_supervised_all:
    input:
        out_files
    output:
        touch(results_dir / "resolvi.done")