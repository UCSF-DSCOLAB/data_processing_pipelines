# What

This directory is meant to contain single cell pre-processing code using the snakemake workflow management system. 

Snakemake is being evaluated as as alternative to nextflow.

## Depooling Pipepline

This pipeline de-pools de-multiplexed `.rds` files that have been de-multiplexed via `freemuxlet` or `demuxlet`.
The input to the pipeline is a set of multiplexed `.rds` files, ie: `TEST-POOL-DM1-SCG1_raw.rds` and the output is a  
number of de-poooled `.rds` files. Each `.rds` file contains a de-pooled object. Files are named by the cluster id or
the bio-specimen id by which they were de-multiplexed.

## Usage

### Step 1 - Module Files

You must first add `/krummellab/data1/modulefiles` to the env var: `MODULEPATH`.

Confirm it works: `module avail`

You should see the mamba module

```
 module avail

------------------------------------------------- /software/c4/modulefiles/repos --------------------------------------------------
   CBI    UserContributed    WitteLab

-------------------------------------------------- /krummellab/data1/modulefiles --------------------------------------------------
   mamba 

```

I've created an environment module called mamba, that points to a mamba installation on c4. Mamba is 
the recommended way to interact and run snakemake.

### Step 2 - Config

Add your de-multiplexed files as input to: `pipelines/depooling/config.yaml`.

### Step 3 - Invoke snakemake

To run a pipeline, use the helper bash script: `./run_snakemake.sh`.

This script:
- Sets up  a designated `mamba` environment that snakemake runs in. 
- Invokes `run_pipeline.py`

Modify the arguments to `run_pipeline.py` as you like. (See below)


### Step 4 - Set arguments for `run_pipeline.py`

This is just a python script that builds a snakemake command.

`./run_pipeline.py --help`

For now this command will look like:

`python run_pipeline.py depooling --singularity-args "--bind /krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/" --local-cores 2`

#### Bind mounts

In one of the commands above you will notice: `"--bind /krummellab/data1/{USER}/tutorial_lib_sep/data/single_cell_GEX/processed/"`.
This mounts the host filesystem into the container filesystem, and note that there is no explicit filesystem target
ie: `"--bind /krummellab/data1/{USER}/tutorial_lib_sep/data/single_cell_GEX/processed/:/container-filesystem/target/"`.
As such the host path gets copied directly into the container path. This also simplifies many input/output specifications
in our snakefile.

### What about slurm?

There is currently a bug in the way snakamake interacts with singularity images via slurm. I've opened an issue for this:
https://github.com/snakemake/snakemake/issues/2808. Until this is resolved, you must invoke this pipeline locally.

## Snakemake

There are a few breaking change from snakemake v7 to v8, most notably (and relevant to us) how we configure snakemake
to interact with slurm. Newer versions of snakemake require the use of plugins to interact with external batch systems:
https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html.

### Custom resource definitions for depooling

If you want to submit custom resource definitions to slurm, you must use:

`python run_pipeline.py depooling --workflow-profile "profiles/pipelines/depooling/"`

You will need to change the config file in the `profiles/pipelines/depooling/`
This will overwrite the standard profile in `profiles/generic`. 

### Tests 

For now we have a very simple (hard coded) script to test de-pooling located at `bin/test_de_pool.R`. 

To test:

- Run the de-pooling pipeline
- Modify the hardcoded paths in `bin/test_de_pool.R`
- Run `bin/test_de_pool.R`