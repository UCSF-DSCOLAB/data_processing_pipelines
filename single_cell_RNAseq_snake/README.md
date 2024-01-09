# What

This directory is meant to contain single cell pre-processing code using the snakemake workflow management system. 

Snakemake is being evaluated as as alternative to nextflow.

## Snakemake set-up

### How was snakemake installed? 

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

### Pipeline execution and mount points

In general, to run a pipeline:

- `conda activate snakemake` to enter the conda snakemake virtual environment
- from `single_cellRNAseq_snake/` , run: `snakemake -s pipelines/de-pooling/Snakefile --cores 1 --use-singularity --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/data/"`

In the command above you will notice: `"--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/data/"`. This mounts the host filesystem into the container filesystem.

An important part of snakemake is the `input` and `output` directives.
  - Snakemake requires files in `input` to be exists when a rule is run, and when the rule finishes, the file(s) in `output` must also exist
  - So you must specify these values in the `config.yaml`

However, since these files are mounted into the container under a different path `/data/{POOL}/automated_processing/file.rds`, 
we cannot use the input/output host file system paths in our R script . Therefore, we have two variables in our `config.yaml`, 
`container_input_files` and `container_output_file` where you specify the container path.

## Depooling

This pipeline de-pools de-multiplexed `.rds` files that have been de-multiplexed via `freemuxlet` or `demuxlet`.
The input to the pipeline is a set of multiplexed `.rds` files, ie: `TEST-POOL-DM1-SCG1_raw.rds` and the output is a 
de-pooled `.rds` file (within contains a list of de-pooled objects).

### Execution

- Add your de-multiplexed files as input to: `pipelines/de-pooling/config.yaml`. (Note see: `example-config/config.yaml` as a reference)
- `conda activate snakemake` to enter the conda snakemake virtual environment.
- from `single_cellRNAseq_snake/` , run: `snakemake -s pipelines/de-pooling/Snakefile --cores 1 --use-singularity --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/data/"`

### Tests (WIP)

This is an active area of research. 

For now we have a simple script to test de-pooling located at `bin/test_de_pool.R`. 

To test:

- Run the de-pooling pipeline
- Modify the hardcoded paths in `bin/test_de_pool.R`
- Run `bin/test_de_pool.R`