from pathlib import Path

# cfg paths
scrnaseq_dir = Path(config['scrnaseq_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
out_files = []

for scrnaseq_reference_rds in (scrnaseq_references := scrnaseq_dir.iterdir()):
    name = scrnaseq_reference_rds.stem

    out_dir = results_dir / f'seurat_to_h5/{name}' 
    out_file = out_dir / '.done'
    out_files.append(out_file)

    rule:
        name: f'seurat_to_h5/{name}'
        input:
            scrnaseq_reference_rds=scrnaseq_reference_rds,
        output:
            out_file=out_file
        params:
            out_dir=out_dir,
        threads: 1
        resources:
            mem='80GB',
            runtime='2h',
        conda:
            "rctd"
        shell:
            """
            Rscript --vanilla workflow/scripts/scRNAseq/seurat_to_h5.R \
            {input.scrnaseq_reference_rds} \
            {params.out_dir} \

            echo "DONE"
            """

rule seurat_to_h5_all:
    input:
        out_files