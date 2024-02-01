import subprocess
import click
import os

DEPOOLING = "depooling"
AVAILABLE_PIPELINES = [DEPOOLING]

@click.command()
@click.argument('pipeline')
@click.option('--local', is_flag=True, default=False, help='Do not use Slurm for job scheduling')
@click.option('--profile', default='profiles/generic/', show_default=True, type=str, help='Full path to the Snakemake profile directory')
@click.option('--workflow-profile', default='', show_default=True, type=str, help='Full path to workflow specific profiles - must be directory.')
@click.option('--singularity-args', default='', type=str, help='Additional arguments to pass to Singularity')
def run_snakemake(pipeline, local, profile, workflow_profile, singularity_args):

    if pipeline not in AVAILABLE_PIPELINES:
        print(f"Please enter a valid pipeline, you've entered: {pipeline}")

    snakemake_cmd = f"snakemake"

    # Use Slurm unless --local is specified
    if local:
        snakemake_cmd += " --cores 1"
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
    # python run_pipeline.py depooling
    # python run_pipeline.py depooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/"
    # python run_pipeline.py depooling --workflow-profile "profiles/pipelines/depooling/"
    run_snakemake()