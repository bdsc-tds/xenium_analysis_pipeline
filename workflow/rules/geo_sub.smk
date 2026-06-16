#######################################
#              Functions              #
#######################################

def get_input2gatherFilesPerSampleForGeoSub(wildcards):
    original_sample_id = get_dict_value(
        config,
        cc.WILDCARDS_NAME,
        cc.WILDCARDS_GEO_SUB_SAMPLES_NAME,
        wildcards.geo_sub_sample_id,
    )
    root_dir = normalise_path(
        f'{config["experiments"][cc.EXPERIMENTS_BASE_PATH_NAME]}/{original_sample_id}',
    )
    versions_file = checkpoints.check10xVersions.get(sample_id=original_sample_id).output[0]
    return {
        "raw_data_dir": root_dir,
        "versions": versions_file,
    }


#######################################
#                Rules                #
#######################################

rule gatherFilesPerSampleForGeoSub:
    input:
        unpack(get_input2gatherFilesPerSampleForGeoSub)
    output:
        temp(f'{config["output_path"]}/geo_sub/geo_sub/_{{geo_sub_sample_id}}.done')
    params:
        prefix=lambda wildcards: wildcards.geo_sub_sample_id,
        out_dir=f'{config["output_path"]}/geo_sub/geo_sub'
    log:
        f'{config["output_path"]}/geo_sub/logs/gatherFilesPerSampleForGeoSub_{{geo_sub_sample_id}}.log'
    resources:
        runtime=60
    shell:
        "python3 workflow/scripts/_data_wrapping/gather_files_for_geo_sub.py "
        "-i {input.raw_data_dir} "
        "-o {params.out_dir} "
        "--prefix {params.prefix} "
        "--versions {input.versions} "
        "-l {log} && "
        "touch {output}"

rule gatherFilesForGeoSub:
    input:
        expand(f'{config["output_path"]}/geo_sub/geo_sub/_{{geo_sub_sample_id}}.done', geo_sub_sample_id=GEO_SUB_SAMPLE_ID)
    output:
        temp(f'{config["output_path"]}/geo_sub/geo_sub/_all_files.txt')
    params:
        out_dir=f'{config["output_path"]}/geo_sub/geo_sub'
    resources:
        runtime=30
    run:
        from pathlib import Path
        out_dir = Path(params.out_dir)
        with open(output[0], 'w', encoding='utf-8') as fh:
            for p in sorted(out_dir.iterdir()):
                if p.is_file() and not p.name.startswith('_') and not p.name.endswith('.done'):
                    fh.write(f'{p.absolute()}\n')

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