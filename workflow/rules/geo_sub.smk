#######################################
#              Functions              #
#######################################

def get_input2gatherFilesPerSampleForGeoSub(wildcards):
    sample_id = extract_layers_from_experiments(
        wildcards.geo_sub_sample_id,
        [0, 1, 2, 3],
        sep_in="_",
        sep_out="/",
        maxsplit=3,
    )[0]
    root_dir = f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{sample_id}'

    return {
        "r_img": f'{root_dir}/morphology.ome.tif',
        "r_ts": f'{root_dir}/transcripts.parquet',
        "p_cnt": f'{root_dir}/cell_feature_matrix.h5',
        "p_cells": f'{root_dir}/cells.parquet',
        "p_cell_boundaries": f'{root_dir}/cell_boundaries.parquet',
        "p_nuc_boundaries": f'{root_dir}/nucleus_boundaries.parquet',
    }


#######################################
#                Rules                #
#######################################

rule gatherFilesPerSampleForGeoSub:
    input:
        unpack(get_input2gatherFilesPerSampleForGeoSub)
    output:
        r_img=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_morphology.ome.tif',
        r_ts=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_transcripts.parquet',
        p_cnt=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cell_feature_matrix.h5',
        p_cells=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cells.parquet',
        p_cell_boundaries=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cell_boundaries.parquet',
        p_nuc_boundaries=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_nucleus_boundaries.parquet'
    log:
        f'{config["output_path"]}/geo_sub/logs/gatherFilesPerSampleForGeoSub_{{geo_sub_sample_id}}.log'
    resources:
        runtime=60
    shell:
        "cp {input.r_img} {output.r_img} && "
        "cp {input.r_ts} {output.r_ts} && "
        "cp {input.p_cnt} {output.p_cnt} && "
        "cp {input.p_cells} {output.p_cells} && "
        "cp {input.p_cell_boundaries} {output.p_cell_boundaries} && "
        "cp {input.p_nuc_boundaries} {output.p_nuc_boundaries} &> {log}"

rule gatherFilesForGeoSub:
    input:
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_morphology.ome.tif', geo_sub_sample_id=GEO_SUB_SAMPLE_ID),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_transcripts.parquet', geo_sub_sample_id=GEO_SUB_SAMPLE_ID),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cell_feature_matrix.h5', geo_sub_sample_id=GEO_SUB_SAMPLE_ID),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cells.parquet', geo_sub_sample_id=GEO_SUB_SAMPLE_ID),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_cell_boundaries.parquet', geo_sub_sample_id=GEO_SUB_SAMPLE_ID),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_sample_id}}_nucleus_boundaries.parquet', geo_sub_sample_id=GEO_SUB_SAMPLE_ID)
    output:
        temp(f'{config["output_path"]}/geo_sub/geo_sub/_all_files.txt')
    resources:
        runtime=30
    run:
        from pathlib import Path
        with open(output[0], 'w', encoding='utf-8') as fh:
            for i in input:
                fh.write(f'{Path(i).absolute()}\n')

rule computeMd5ForGeoSub:
    input:
        f'{config["output_path"]}/geo_sub/geo_sub/_all_files.txt'
    output:
        f'{config["output_path"]}/geo_sub/geo_sub/MD5.txt'
    log:
        f'{config["output_path"]}/geo_sub/logs/computeMd5ForGeoSub.log'
    threads:
        4
    resources:
        mem_mb=lambda wildcards, input, attempt: min(2048 * attempt * 2, 40960)
    shell:
        "mamba run -n utilnest python workflow/scripts/_data_wrapping/compute_checksum.py "
        "--batch_file {input} "
        "-o {output} "
        "--algo md5 "
        "-t {threads} "
        "-l {log}"
