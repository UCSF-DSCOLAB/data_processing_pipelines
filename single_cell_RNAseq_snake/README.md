# What

This directory is meant to contain single cell pre-processing code using the snakemake workflow management system. 

Snakemake it being evaluated as as alternative to nextflow.

## Install and usage

### Miniforge

`mamba` is the recommended way to install snakemake. It is a conda distribution that some optimized features suitable for snakemake. 

Miniforge is a repository that holds the minimal install for Conda.

The conda distribution is installed in `/krummellab/data1/software/miniforge3` and the utilities like `mamba` live in: `/krummellab/data1/software/miniforge3/bin`.

### Usage

Add `/krummellab/data1/software/miniforge3/bin` to your path.

Run `conda activate snakemake` to enter the conda snakemake virtual environment.

Run `snakemake -h` to verify that you can reference snakemake globally.

### What is going on?

- `conda` is a package manager that allows you to manage your dependencies within an isolated environment, very similar to `virtualenv`

- We used `conda` to create a snakemake virtual environment called `snakemake`. This virtual environment lives in: `/krummellab/data1/software/miniforge3/envs/snakemake`, and is essentially an installation of snakemake.

- You can "activate" or rather point to an installation / or conda environment by running `conda activate snakemake`. Once you do so, you can reference `snakemake` as if it were a binary in your path.


## Pipeline execution 

- `conda activate snakemake` to enter the conda snakemake virtual environment
- from `single_cellRNAseq_snake/` , run: `snakemake -s pipelines/de-pooling/Snakefile --cores 1 --use-singularity --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/data/"`

### How do mount points work with snakmake/singularity?

- In the command above you will notice: `"--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/data/"`. This mounts the host filesystem into the container filesystem.
- An important part of snakemake is the `input` and `output` directives.
  - Snakemake requires files in `input` to be exists when a rule is run, and when the rule finishes, the file(s) in `output` must also exist
  - So you must specify these values in the `config.yaml`
- However, since these files are mounted into the container under a different path `/data/{POOL}/automated_processing/file.rds`, we cannot use the input/output host file system paths in our R script
- Therefore, we have two variables in our `config.yaml`, `container_input_files` and `container_output_file` where you specify the container path

## Tests

#TODO


