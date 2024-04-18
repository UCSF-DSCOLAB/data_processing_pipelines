import subprocess
import os
import argparse

DEPOOLING = "depooling"
AVAILABLE_PIPELINES = [DEPOOLING]

def run_snakemake(pipeline, local_cores, profile, workflow_profile, singularity_args):

    if pipeline not in AVAILABLE_PIPELINES:
        print(f"Please enter a valid pipeline, you've entered: {pipeline}")
        return

    snakemake_cmd = "snakemake"

    # Use Slurm unless --local is specified
    if local_cores:
        snakemake_cmd += f" --cores {local_cores}"
    else:
        snakemake_cmd += f" --profile {os.path.abspath(profile)}"
        if workflow_profile:
            snakemake_cmd += f" --workflow-profile {os.path.abspath(workflow_profile)}"

    snakemake_cmd += f" --snakefile pipelines/{pipeline}/Snakefile"

    # Pipeline specific configs
    if pipeline == DEPOOLING:
        # Pass in custom args
        if singularity_args:
            snakemake_cmd += f" --use-singularity --singularity-args '{singularity_args}'"
        else:
            bind_mount_path = f"/krummellab/data1/{os.getenv('USER')}/tutorial_lib_sep/data/single_cell_GEX/processed/"
            snakemake_cmd += f" --use-singularity --singularity-args '--bind {bind_mount_path}'"

    print(f"Running command: {snakemake_cmd}")
    subprocess.run(snakemake_cmd, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Snakemake pipelines.")
    parser.add_argument('pipeline', type=str, help="Pipeline to run")
    parser.add_argument('--local-cores', type=int, default=None, help="Do not use Slurm for job scheduling. Include the number of local cores")
    parser.add_argument('--profile', default='profiles/generic/', type=str, help="Full path to the Snakemake profile directory")
    parser.add_argument('--workflow-profile', default='', type=str, help="Full path to workflow specific profiles - must be directory.")
    parser.add_argument('--singularity-args', default='', type=str, help="Additional arguments to pass to Singularity")

    args = parser.parse_args()
    run_snakemake(args.pipeline, args.local_cores, args.profile, args.workflow_profile, args.singularity_args)
