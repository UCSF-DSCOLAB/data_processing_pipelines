import subprocess
import click
import os

@click.command()
@click.argument('pipeline')
@click.option('--local', is_flag=True, default=False, help='Do not use Slurm for job scheduling')
@click.option('--profile-path', default='slurm/', show_default=True, type=str, help='Full path to the Snakemake profile directory')
@click.option('--singularity-args', default='', type=str, help='Additional arguments to pass to Singularity')
@click.option('--cluster-config', default='pipelines/de-pooling/cluster_config.yaml', type=str, help='Path to the cluster config')
def run_snakemake(pipeline, local, profile_path, singularity_args, cluster_config):

    base_dir = "./"
    snakemake_cmd = f"snakemake -s pipelines/{pipeline}/Snakefile"

    if singularity_args:
        snakemake_cmd += f" --use-singularity --singularity-args '{singularity_args}'"

    # Use Slurm unless --local is specified
    if local:
        snakemake_cmd += " --cores 1"
    else:
        snakemake_cmd += f" --profile {os.path.abspath(profile_path)} --cluster-config {cluster_config}"

    print(f"Running command: {snakemake_cmd}")
    subprocess.run(snakemake_cmd, shell=True)

if __name__ == "__main__":
    # python run_pipeline.py de-pooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/"
    run_snakemake()