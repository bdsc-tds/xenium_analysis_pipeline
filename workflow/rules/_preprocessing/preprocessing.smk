# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# params from existing config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
n_comps = 50
n_neighbors = 50
min_dist = 0.3
metric = 'cosine'

out_files_panel = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):

            k = (segmentation.stem,cohort.stem,panel.stem)
            name = '/'.join(k)

            if replicate_transcripts_path.exists():

                out_file = results_dir / f'embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                out_files_panel.append(out_file)

            rule:
                name: f'embed_panel/{name}'
                input:
                    panel=panel,
                output:
                    out_file=out_file,
                params:
                    n_comps=n_comps,
                    n_neighbors=n_neighbors,
                    metric=metric,
                    min_dist=min_dist,
                    min_counts=min_counts,
                    min_features=min_features,
                    max_counts=max_counts,
                    max_features=max_features,
                    min_cells=min_cells,
                threads: 1
                resources:
                    mem='100GB' if panel.stem == '5k' else '50GB',
                    runtime='30m' if panel.stem == '5k' else '20m',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/_preprocessing/embed_panel.py \
                        --panel {input.panel} \
                        --out_file {output.out_file} \
                        --n_comps {params.n_comps} \
                        --n_neighbors {params.n_neighbors} \
                        --metric {params.metric} \
                        --min_dist {params.min_dist} \
                        --min_counts {params.min_counts} \
                        --min_features {params.min_features} \
                        --max_counts {params.max_counts} \
                        --max_features {params.max_features} \
                        --min_cells {params.min_cells}

                    echo "DONE"
                    """


out_files_cohort = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 

            k = (segmentation.stem,cohort.stem)
            name = '/'.join(k)

            if replicate_transcripts_path.exists():

                out_file = results_dir / f'embed_cohort/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                out_files_cohort.append(out_file)

            rule:
                name: f'embed_cohort/{name}'
                input:
                    cohort=cohort,
                output:
                    out_file=out_file,
                params:
                    n_comps=n_comps,
                    n_neighbors=n_neighbors,
                    metric=metric,
                    min_dist=min_dist,
                    min_counts=min_counts,
                    min_features=min_features,
                    max_counts=max_counts,
                    max_features=max_features,
                    min_cells=min_cells,
                threads: 1
                resources:
                    mem='100GB' if panel.stem == '5k' else '50GB',
                    runtime='30m' if panel.stem == '5k' else '20m',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1'
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/_preprocessing/embed_cohort.py \
                        --cohort {input.cohort} \
                        --out_file {output.out_file} \
                        --n_comps {params.n_comps} \
                        --n_neighbors {params.n_neighbors} \
                        --metric {params.metric} \
                        --min_dist {params.min_dist} \
                        --min_counts {params.min_counts} \
                        --min_features {params.min_features} \
                        --max_counts {params.max_counts} \
                        --max_features {params.max_features} \
                        --min_cells {params.min_cells}

                    echo "DONE"
                    """

rule embed_panel_samples:
    input:
        out_files_panel

rule embed_cohort_samples:
    input:
        out_files_cohort