#!/bin/bash
#SBATCH --job-name=bulk_rnaseq_smk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=/krummellab/data1/%u/logs/bulk_rnaseq_smk_%j.log
#SBATCH --exclude=c4-n20

# Usage:
#   sbatch ./run_snakemake_c4_singularity.sh
# Forward extra Snakemake args to this script after "--":
#   sbatch ./run_snakemake_c4_singularity.sh -- --rerun-incomplete -k

set -euo pipefail

# Singularity image (same as Nextflow container.config)
SIF="/krummellab/data1/singularity_images/bulk_rnaseq/v6/bulk_rna_seq.sif"

# Scratch working directory (cleanup on success)
work_base="/c4/scratch/${USER}/snakemake"
work="${work_base}/${SLURM_JOB_ID}"
mkdir -p "${work}"

# Propagate tmpdir into the container
export TMPDIR="${work}"
export SINGULARITYENV_TMPDIR="${work}"
export APPTAINERENV_TMPDIR="${work}"

# Cache for Singularity to avoid filling $HOME
export SINGULARITY_CACHEDIR="${work_base}/cache"
export APPTAINER_CACHEDIR="${work_base}/cache"

# Bind mounts (adjust as needed)
BIND_DIRS="$(pwd),/krummellab,/c4,/scratch,/home/${USER}"

# Respect Slurm CPUs per task for Snakemake --cores
CORES="${SLURM_CPUS_PER_TASK:-16}"

cleanup() {
  status=$?
  if [[ $status -eq 0 ]]; then
    echo "Pipeline complete; cleaning tmp ${work}"
    rm -rf "${work}"
  else
    echo "Pipeline failed; keeping tmp at ${work}"
  fi
}
trap cleanup EXIT

# Optional: verify Snakemake is available inside the container
if ! singularity exec --cleanenv -B "${BIND_DIRS}" "${SIF}" bash -lc 'command -v snakemake >/dev/null 2>&1'; then
  echo "ERROR: snakemake is not available inside ${SIF}."
  echo "Either install Snakemake into the image or use the conda-based script run_snakemake_c4.sh."
  exit 1
fi

# Run Snakemake inside the container. No --use-conda since tools should be in the SIF.
singularity exec --cleanenv -B "${BIND_DIRS}" "${SIF}" \
  snakemake -s Snakefile \
    --cores "${CORES}" \
    --latency-wait 120 \
    --printshellcmds \
    --directory "$(pwd)" \
    "$@"