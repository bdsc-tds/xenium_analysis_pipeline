import os
import glob
import pandas as pd
import config_constants as cc

# optional mapping CSV
mapping_file = config.get('name_mapping_csv')
brcode_to_ptcode = {}

if mapping_file and os.path.exists(mapping_file):
    try:
        df_mapping = pd.read_csv(mapping_file, dtype=str)
        brcode_to_ptcode = dict(zip(df_mapping.iloc[:,0], df_mapping.iloc[:,1]))
    except Exception as e:
        print(f"Warning: Could not load mapping file {mapping_file}: {e}")

# longest-substring matching
# e.g. ensures "0WJ3_big" is replaced before "0WJ3"
sorted_brcodes = sorted(brcode_to_ptcode.keys(), key=len, reverse=True)

# map renamed flattened IDs -> original paths
base_path = config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]
search_pattern = os.path.join(base_path, "*", "*", "*", "*") # condition/panel/donor/sample
found_dirs = glob.glob(search_pattern)

target_id_to_path = {}

for d in found_dirs:
    rel_path = os.path.relpath(d, base_path).replace("\\", "/")
    final_id = rel_path.replace("/", "_")
    
    for brcode in sorted_brcodes:
        if brcode in final_id:
            final_id = final_id.replace(brcode, brcode_to_ptcode[brcode])
            break

    target_id_to_path[final_id] = rel_path

# final list of wildcards
TARGET_GEO_IDS = list(target_id_to_path.keys())


#######################################
#              Functions              #
#######################################

def get_input2gatherFilesPerSampleForGeoSub(wildcards):
    original_path = target_id_to_path[wildcards.geo_sub_id]
    
    root_dir = f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{original_path}'

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
        r_img=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_morphology.ome.tif',
        r_ts=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_transcripts.parquet',
        p_cnt=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cell_feature_matrix.h5',
        p_cells=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cells.parquet',
        p_cell_boundaries=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cell_boundaries.parquet',
        p_nuc_boundaries=f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_nucleus_boundaries.parquet'
    log:
        f'{config["output_path"]}/geo_sub/logs/gatherFilesPerSampleForGeoSub_{{geo_sub_id}}.log'
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
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_morphology.ome.tif', geo_sub_id=TARGET_GEO_IDS),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_transcripts.parquet', geo_sub_id=TARGET_GEO_IDS),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cell_feature_matrix.h5', geo_sub_id=TARGET_GEO_IDS),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cells.parquet', geo_sub_id=TARGET_GEO_IDS),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_cell_boundaries.parquet', geo_sub_id=TARGET_GEO_IDS),
        expand(f'{config["output_path"]}/geo_sub/geo_sub/{{geo_sub_id}}_nucleus_boundaries.parquet', geo_sub_id=TARGET_GEO_IDS)
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
    conda:
        "../envs/geo_sub.yml"
    shell:
        "python3 workflow/scripts/_data_wrapping/compute_checksum.py "
        "--batch_file {input} "
        "-o {output} "
        "--algo md5 "
        "-t {threads} "
        "-l {log}"