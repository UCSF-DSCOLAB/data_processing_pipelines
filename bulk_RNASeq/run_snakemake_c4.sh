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
#   sbatch ./run_snakemake_c4.sh
# Forward extra Snakemake args to this script after "--":
#   sbatch ./run_snakemake_c4.sh -- --rerun-incomplete -k

set -euo pipefail

# Singularity image (same as Nextflow container.config)
SIF="/krummellab/data1/singularity_images/bulk_rnaseq/v6/bulk_rna_seq.sif"

# Scratch working directory (cleanup on success)
work_base="/c4/scratch/${USER}/snakemake"
work="${work_base}/${SLURM_JOB_ID}"
mkdir -p "${work}"

# Use a scratch HOME inside the container to avoid failing mounts
HOME_DIR="${work_base}/home"
mkdir -p "${HOME_DIR}"

# Propagate tmpdir into the container (Apptainer preferred; keep Singularity for compatibility)
export TMPDIR="${work}"
export APPTAINERENV_TMP="${work}"
export APPTAINERENV_TMPDIR="${work}"
export SINGULARITYENV_TMP="${work}"
export SINGULARITYENV_TMPDIR="${work}"
# Avoid leaking host user site-packages into container
export APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_PYTHONNOUSERSITE=1

# Cache for Singularity to avoid filling $HOME
export SINGULARITY_CACHEDIR="${work_base}/cache"
export APPTAINER_CACHEDIR="${work_base}/cache"
mkdir -p "${SINGULARITY_CACHEDIR}"

# Bind mounts
HOST_PWD="$(pwd)"
BIND_DIRS="${HOST_PWD}:/work,/krummellab,/c4,/scratch"
# Only bind host home if it actually exists on this node
if [[ -d "/c4/home/${USER}" ]]; then
  BIND_DIRS="${BIND_DIRS},/c4/home/${USER}"
fi

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

# Optional: verify container launches and snakemake is available
if ! singularity exec --cleanenv --home "${HOME_DIR}" -B "${BIND_DIRS}" --pwd /work "${SIF}" \
     bash -lc 'command -v snakemake >/dev/null 2>&1'; then
  echo "ERROR: snakemake is not available inside ${SIF} (or container failed to start)."
  exit 1
fi

# Run Snakemake inside the container. No --use-conda since tools should be in the SIF.
singularity exec --cleanenv --home "${HOME_DIR}" -B "${BIND_DIRS}" --pwd /work "${SIF}" \
  snakemake -s Snakefile \
    --cores "${CORES}" \
    --latency-wait 120 \
    --printshellcmds \
    --directory /work \
    "$@"
