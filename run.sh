#!/bin/bash


##############################
#         USER SETUP         #
##############################

# For slurm clusters: modules to be loaded
MODULES="singularityce"

# The name of or path to either `mamba` or `conda`.
CONDA_BIN=mamba

# The name of or path to the conda environment.
ENV_NAME=xenium_analysis_pipeline

# Path to the profile for submitting jobs under "local" mode.
LOCAL_PROFILE=profiles/local

# Path to the profile for submitting jobs under "cluster" mode.
CLUSTER_PROFILE=profiles/slurm

# Array of directories to bind to containers.
# Each element should be in the following form: LOCAL_DIR:SINGULARITY_DIR
# Inexistent local directories will be filtered out.
SINGULARITY_BIND_DIRS=(  )

# Temporary folder for Snakemake.
# Leave blank to use the default location.
SNAKEMAKE_CACHE_DIR=


###############################
#         DO NOT EDIT         #
###############################

# Help message.
help()
{
    echo "Usage: [ -m | --mode MODE ] [ -c | --core CORE ] [ -n | --dry-run ] [ -R | --forcerun RULE ] [ -U | --until RULE ] [ --dag OUTPUT ] [ --unlock ] [ -v | --verbose ] [ -h | --help ]
        -m,--mode MODE: the pipeline will be run on 'local' (default) or on 'cluster'.
        -c,--core CORE: the number of cores to be used when -m,--mode is unset or 'local' (default: 1); ignored when -m,--mode is 'cluster'.
        -j,--jobs JOBS: the number of jobs submitted to the cluster at the same time when -m,--mode is 'cluster'. (default: 500).
        --retries RETRIES: the number of retries for failed jobs. (default: 0).
        -n,--dry-run: dry run.
        -R,--forcerun RULE: force the re-execution or creation of the given rule or file. Repeat this option multile times for multiple rules or files.
        -U,--until RULE: runs the pipeline until it finishes the specified rule or generated the file. Repeat this option multile times for multiple rules or files.
        --dag OUTPUT: draw dag and save to OUTPUT.pdf.
        --unlock: unlock the working directory.
        -v,--verbose: print more information.
        -h,--help: print this message.
        "
    
    exit "$1"
}

# Allowed command line arguments.
SHORT_OPTS=m:,c:,j:,n,R:,-U:,v,h
LONG_OPTS=mode:,core:,jobs:,retries:,dry-run,forcerun:,until:,dag:,unlock,verbose,help
OPTS=$(getopt -n xenium_analysis_pipeline --options $SHORT_OPTS --longoptions $LONG_OPTS -- "$@")
eval set -- "$OPTS"

# Bind directories, if they are present.
IFS=':'
for i in "${SINGULARITY_BIND_DIRS[@]}"
do
    :
    read -r -a paths <<< "$i"
    if test -x "${paths[0]}"; then
        if [ -z "${SINGULARITY_BIND}" ]; then
            SINGULARITY_BIND="${i}"
        else
            SINGULARITY_BIND="${SINGULARITY_BIND},${i}"
        fi
    fi
done
IFS=''

# Constants.
SINGULARITY_BIND_OPT="--bind $SINGULARITY_BIND"
OTHER_OPT=(--rerun-triggers mtime --rerun-incomplete --software-deployment-method conda apptainer --apptainer-args "--nv --no-home --cleanenv --env RUST_BACKTRACE=full $SINGULARITY_BIND_OPT" -kp)

# Options related to $CONDA_BIN
if [[ "$CONDA_BIN" = mamba ]]; then
    CONDA_OPT=
elif [[ "$CONDA_BIN" = conda ]]; then
    CONDA_OPT="--no-capture-output"
fi

# Variables.
LOCAL=1
EXEC_OPT=(--profile "$LOCAL_PROFILE")
CORE_OPT=(--cores 1)
JOBS_OPT=(--jobs 500)
RETRIES_OPT=(--retries 0)
DRY_RUN_OPT=
FORCE_RUN_OPT=()
UNTIL_OPT=()
DAG_OPT=
UNLOCK=0
VERBOSE=0

# Process command line arguments.
while :
do
    case "$1" in

        -m | --mode)
            case "$2" in

                local)
                    ;;

                cluster)
                    LOCAL=0
                    EXEC_OPT=(--profile "$CLUSTER_PROFILE")
                    ;;

                *)
                    echo "Invalid value for -m,--mode: $2"
                    help 1
                    ;;

            esac
            shift 2
            ;;

        -c | --core)
            if [[ $2 =~ ^[1-9][0-9]*$ ]]; then
                CORE_OPT=(--cores "$2")
            else
                echo "Invalid value for -c,--core: $2"
                help 1
            fi
            shift 2
            ;;

        -j | --jobs)
            if [[ $2 =~ ^[1-9][0-9]*$ ]]; then
                JOBS_OPT=(--jobs "$2")
            else
                echo "Invalid value for -j,--jobs: $2"
                help 1
            fi
            shift 2
            ;;
        
        --retries)
            if [[ $2 =~ ^[0-9]+$ ]]; then
                RETRIES_OPT=(--retries "$2")
            else
                echo "Invalid value for --retries: $2"
                help 1
            fi
            shift 2
            ;;

        -n | --dry-run)
            DRY_RUN_OPT=-n
            shift 1
            ;;

        -R | --forcerun)
            if [ -z "$2" ]; then
                echo "Empty value for -R,--forcerun."
                help 1
            fi

            FORCE_RUN_OPT=("${FORCE_RUN_OPT[@]}" "$2")

            shift 2
            ;;

        -U | --until)
            if [ -z "$2" ]; then
                echo "Empty value for -U,--until."
                help 1
            fi

            UNTIL_OPT=("${UNTIL_OPT[@]}" "$2")

            shift 2
            ;;

        --dag)
            DAG_OPT="$2"

            if [[ -z $DAG_OPT ]]; then
                echo "Empty value for --dag."
                help 1
            fi

            break
            ;;

        --unlock)
            UNLOCK=1
            break
            ;;
        
        -v | --verbose)
            VERBOSE=1
            shift 1
            ;;

        -h | --help)
            help 0
            ;;

        --)
            shift;
            break
            ;;

        *)
            echo "Unknown option: $1"
            help 1
            ;;

    esac
done

# Set Snakemake runtime temporary directory
if [[ -n "$SNAKEMAKE_CACHE_DIR" ]]; then
    export XDG_CACHE_HOME="$SNAKEMAKE_CACHE_DIR"
fi

# Draw dag and save to disk.
# Priority: 1 (highest; other options will be ommitted)
if [[ -v DAG_OPT && -n $DAG_OPT ]]; then
    $CONDA_BIN run $CONDA_OPT -n $ENV_NAME snakemake --dag | dot -Tpdf > $DAG_OPT.pdf
    exit 0
fi

# Unlock the working directory.
# Priority: 2 (2nd highest; other options will be ommitted)
if [[ $UNLOCK -eq 1 ]]; then
    $CONDA_BIN run $CONDA_OPT -n $ENV_NAME snakemake --unlock
    exit 0
fi

# Command for snakemake.
COMPLETE_CMD=(snakemake "${OTHER_OPT[@]}")

# Verbose for snakemake.
if [[ $VERBOSE -eq 1 ]]; then
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" "--verbose")
fi

# Profile for job submission.
COMPLETE_CMD=("${COMPLETE_CMD[@]}" "${EXEC_OPT[@]}")

# Postprocess after parsing command line arguments.
if [[ $LOCAL -eq 0 ]]; then
    # On cluster.
    CORE_OPT=()

    # Number of jobs to be submitted to the cluster.
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" "${JOBS_OPT[@]}")
else
    # On local.
    JOBS_OPT=()

    # Number of local cores to be used for snakemake.
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" "${CORE_OPT[@]}")
fi

# Dry run mode.
# Priority: 3
if [[ -n "$DRY_RUN_OPT" ]]; then
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" "$DRY_RUN_OPT")
fi

# Force run rules.
# Priority: 3
if [[ -n "${FORCE_RUN_OPT[*]}" ]]; then
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" --forcerun "${FORCE_RUN_OPT[@]}")
fi

# Stop after rules.
# Priority: 3
if [[ -n "${UNTIL_OPT[*]}" ]]; then
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" --until "${UNTIL_OPT[@]}")
fi

# Retries for failed jobs (not for dry-run mode).
if [[ -z "$DRY_RUN_OPT" ]]; then
    COMPLETE_CMD=("${COMPLETE_CMD[@]}" "${RETRIES_OPT[@]}")
fi

if [[ $VERBOSE -eq 1 ]]; then
    printf "> Command used to run snakemake:\n  >> "
    echo "${COMPLETE_CMD[@]}"
fi

# Load packages on cluster.
module load $MODULES

# Run snakemake command along with options.
$CONDA_BIN run $CONDA_OPT -n $ENV_NAME "${COMPLETE_CMD[@]}"
