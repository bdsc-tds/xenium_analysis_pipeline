from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
methods = ['conditional','jaccard','pearson','spearman']
target_counts = [30,50,200]
out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):
            for sample in (samples := panel.iterdir()):
                for replicate in (replicates := sample.iterdir()):

                    k = (segmentation.stem,cohort.stem,panel.stem,sample.stem,replicate.stem)
                    replicate_path = replicate / "normalised_results/outs"
                    name = '/'.join(k)

                    if replicate_path.exists():

                        for method in methods:
                            for target_count in target_counts:
                                if target_count > 50 and panel.stem != '5k':
                                    continue

                                out_file_coexpr = results_dir / f'coexpression/{name}/coexpression_{method}_{target_count}.parquet' 
                                out_file_pos_rate = results_dir / f'coexpression/{name}/positivity_rate_{method}_{target_count}.parquet'

                                out_files.extend([out_file_coexpr,out_file_pos_rate])

                                rule:
                                    name: f'coexpression/{name}/{method}_{target_count}'
                                    input:
                                        replicate_path=replicate_path,
                                    output:
                                        out_file_coexpr=out_file_coexpr,
                                        out_file_pos_rate=out_file_pos_rate,
                                    params:
                                        method=method,
                                        target_count=target_count,
                                    threads: 1
                                    resources:
                                        mem='40GB' if panel.stem == '5k' else '20GB',
                                        runtime='20m' if panel.stem == '5k' else '10m',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file_coexpr})"

                                        python workflow/scripts/xenium/coexpression_sample.py \
                                        {input.replicate_path} \
                                        {output.out_file_coexpr} \
                                        {output.out_file_pos_rate} \
                                        {params.method} \
                                        {params.target_count} \

                                        echo "DONE"
                                        """


rule coexpression_samples:
    input:
        out_files