#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --exclude=c4-n20

# to run:
#   sbatch ./run_pipeline_c4.sh
# to resume:
#.  sbatch ./run_pipeline_c4.sh -resume

# Arugments:
# pass as many additional arguments to nextflow as you'd like (e.g. -resume, -with-timeline, -profile test)

nextflow_args=("$@")
results_directory=""

# A command-line override takes precedence over config/user_directories.config.
for ((i = 0; i < ${#nextflow_args[@]}; i++)); do
    case "${nextflow_args[$i]}" in
        --results_directory)
            results_directory="${nextflow_args[$((i + 1))]:-}"
            ;;
        --results_directory=*)
            results_directory="${nextflow_args[$i]#*=}"
            ;;
    esac
done

if [ -z "${results_directory}" ]; then
    results_directory=$(sed -nE 's/^[[:space:]]*results_directory[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p' config/user_directories.config)
    results_directory=${results_directory//'${USER}'/${USER}}
fi

if [[ -z "${results_directory}" || "${results_directory}" != /* ]]; then
    echo "results_directory must be configured as an absolute path" >&2
    exit 1
fi

mkdir -p "${results_directory}/logs"
exec > >(tee -a "${results_directory}/logs/bulk_rnaseq_nf_${SLURM_JOB_ID}.log") 2>&1

# check on arguments
failed=false


if [ "$failed" = true ];
then
  exit 1
else
    # create a working directory
    nf_work=${results_directory}/nextflow/
    mkdir -p "${nf_work}"
    export NXF_WORK="${nf_work}"
    export APPTAINERENV_TMPDIR="${nf_work}"
    # run the pipeline
    nextflow run bulk_rna_seq.nf -c config/base.config -profile hpc_c4 "${nextflow_args[@]}"
fi
