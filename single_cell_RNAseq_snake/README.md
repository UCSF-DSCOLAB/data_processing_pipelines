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

Note: this conda environment uses snakemake `8.4.2`

### What is going on?

- `conda` is a package manager that allows you to manage your dependencies within an isolated environment, very similar to `virtualenv`

- We used `conda` to create a snakemake virtual environment called `snakemake`. This virtual environment lives in: `/krummellab/data1/software/miniforge3/envs/snakemake`, and is essentially an installation of snakemake.

- You can "activate" or rather point to an installation / or conda environment by running `conda activate snakemake`. Once you do so, you can reference `snakemake` as if it were a binary in your path.

### Pipeline execution and mount points

To run a pipeline, use the helper python script: `run_pipeline.py`. For usage: `python run_pipeline.py --help`.

In general:
- `conda activate snakemake` to enter the conda snakemake virtual environment 
- `python run_pipeline.py depooling`


### Snakemake 8.4.2

There are a few breaking change from snakemake v7 to v8, most notably (and relevant to us) how we configure snakemake
to interact with slurm. Newer versions of snakemake require the use of plugins to interact with external batch systems:
https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html.

## Depooling

This pipeline de-pools de-multiplexed `.rds` files that have been de-multiplexed via `freemuxlet` or `demuxlet`.
The input to the pipeline is a set of multiplexed `.rds` files, ie: `TEST-POOL-DM1-SCG1_raw.rds` and the output is a  
number of de-poooled `.rds` files. Each `.rds` file contains a de-pooled object. Files are named by the cluster id or
the bio-specimen id by which they were de-multiplexed.

### Execution

- Add your de-multiplexed files as input to: `pipelines/depooling/config.yaml`. 
- `conda activate snakemake` to enter the conda snakemake virtual environment.
- from `single_cellRNAseq_snake/`, `python run_pipeline depooling`

Some other notable parameters:

- `python run_pipeline.py depooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/"`
- `python run_pipeline.py depooling --workflow-profile "profiles/pipelines/depooling/"`

In one of the commands above you will notice: `"--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/"`.
This mounts the host filesystem into the container filesystem, and note that there is no explicit filesystem target
ie: `"--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/:/container-filesystem/target/"`.
As such the host path gets copied directly into the container path. This also simplifies many input/output specifications
in our snakefile.

### Custom resource definitions for depooling

If you want to submit custom resource definitions to slurm, you must use:

`python run_pipeline.py depooling --workflow-profile "profiles/pipelines/depooling/"`

You will need to change the config file in the `profiles/pipelines/depooling/`
This will overwrite the standard profile in `profiles/generic`. 

### Tests 

TODO: use snakemake to generate unit tests.

For now we have a very simple (hard coded) script to test de-pooling located at `bin/test_de_pool.R`. 

To test:

- Run the de-pooling pipeline
- Modify the hardcoded paths in `bin/test_de_pool.R`
- Run `bin/test_de_pool.R`