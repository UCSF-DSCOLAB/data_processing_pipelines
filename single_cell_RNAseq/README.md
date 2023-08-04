### Update: 07/17/2023

This code was initially developed at https://github.com/UCSF-DSCOLAB/sc_seq_nextflow/. Please note that the original
repository is no longer maintained and has been marked as deprecated. Going forward, all enhancements, bug fixes, 
and active development are to be carried out in this repository.

## How to

### Running on c4

#### Pre-reqs

In order to run nextflow we need to install a compatible version of java. We unfortunately can't use a package manager
(because of sudo restrictions), so we must download it manually.

1. In your home directory: `mkdir jdk`
2. `wget https://download.oracle.com/java/17/archive/jdk-17.0.7_linux-x64_bin.tar.gz`
3. `mv jdk-17.0.7_linux-x64_bin.tar.gz jdk/`
4. `cd jdk`
5. `tar -xf jdk-17.0.7_linux-x64_bin.tar.gz`
6. Add the following to your `~/.bashrc`: 

```bash
# global environment variables
export JAVA_HOME="/c4/home/${USER}/jdk/jdk-17.0.6"
export NXF_JAVA_HOME="/c4/home/${USER}/jdk/jdk-17.0.6"
```

Install nextflow

1. In your home directory: `mkdir bin` (If it doesn't exist)
2. `cd bin`
3. `wget -qO- https://get.nextflow.io | bash`
4. `chmod +x nextflow`
5. Make sure `${USER}/bin` is added to your path

#### What do I need to configure?

Unfortunately there is a ton of input that is required for pipelines to run. Such as:

- Directories with fastq files
- Location of reference genomes
- Process specific settings and flags

In order to view what all of these settings are, you can check out `nextflow.config`. 
To actually supply the parameters to this pipeline, you must submit a json file with these values.
Some examples include: `example-inputs/param_1.json` `example-inputs/param_2_v2.json`.

There is also a directory called `config/` however these are settings that specific to c4 and typically
do not need to be tweaked. `nextflow.config` imports these settings for you.

#### Data Generation

For those that are not familiar with sequencing and data generation, recall that:

1. We obtain $s$ biological samples / biospecimens. Each sample maps to a unique individual and  
contains $c$ cells.
2. For each sample we isolate $n$ cells, and then load those cells into the sequencer. When $s > 1$ this is
called pooling. 
3. The sequencer runs and generates a library. Each library can contain $r$ reads.
4. We repeat steps 2-3, except we load in a new set of cells.

Note: Pooling is very specific to Colabs

#### What is the input and output of this pipeline?

TODO -> Complete once the input parameters are refactored

## Conventions


### Params

If you have a process that requires input from the `params` keyword, do NOT parse/extract it
at the beginning of the pipeline and supply those inputs to processes downstream. 

Ex:
```
some_param = params.collectMany{ ... }
another_param = params.collectMany { ... }

PROCESS_1()

PROCESS_2(another_param)

PROCESS_3(some_param)
```

Do this instead:

```
// Extract library directories and their corresponding data types
some_param = params.collectMany { ... }

PROCESS_1(some_param)

another_param = params.collectMany { ... }

PROCESS_2(another_param)
```

It is much easier to determine which params belong to which processes this way.

If your params require heavy parsing and manipulation it is probably a sign your 
param structure must be re-designed. 

### Variable names

Groovy is JVM-esque so it might have made sense to use `camelCase` for variable names.
However, most of the code written in this repo is `snake_case`, so lets continue with this
convention.

It is also helpful to prefix channels with the `ch_` prefix.


### Initial SC-Seq Pipeline in Nextflow

Update 1/25/2023: Separated the pipeline into two steps for pre+post QC.


This is a nextflow implementation of the `sicca_sc-seq_pipeline` written by Kim Taylor and Ravi Patel, with the following modifications: (1) addition of a freemuxlet merge step to combine across pools, (2) addition of an option to filter based on QC pre-freemuxlet, and (3) making inclusion of VDJ data optional. Future efforts will make this more flexible for different inputs.

Note: this pipeline assumes the project structure we use in `immunox`

To run:
 `sbatch run.sh <path_to_config.json> <step>`

To resume a previous run:
 `sbatch run_resume.sh <path_to_config.json> <step> <job_id_for_original_run>`

The `step` must be one of ['pre_qc', 'post_qc', pre_fmx_qc', 'post_fmx_qc'].
The pipeline is designed to be run either in the standard way, with (1) pre_qc, (2) setting cutoffs, then (3) post_qc.
If filtering before freemuxlet is desired, such as for tissue data, instead do: (1) pre_fmx_qc, (2) setting cutoffs, then (3) post_fmx_qc.

See `nextflow_schema.json` or `example/library.json` for how to set up the config files.

By default, the nextflow working directory is:
`/c4/scratch/<user>/nextflow/<original_job_id>`
This is deleted on successful completion of an initial run. If the pipeline fails, you can resume it with `run_repeat.sh` (which uses the same working directory), or you should manually remove this directory. After `run_repeat.sh`, you must manually remove this directory. 

#### Pre-QC Pipeline

Please note that the `pipeline_pre_qc.nf` uses a different json structure. 
You can see an example in `example-inputs/param_2_v2.json`. Eventually all pipelines will use
this structure.

#### Next steps:
Below are outlined the next steps for turning this into a usable pipeline

key additions:

- input validation
- improve toy/test cases
- add additional outputs to final seurat step
- create a new singularity container with all R packages needed (e.g. doubletfinder, dsb)
- dynamic memory allotment, retries


missing steps:

- add in freemuxlet assign to genotype step
- add demuxlet step
- ?separate ADT step from seurat_qc?
- ?Harmony/RPCA/WNN?

parameters:

- add options to keep.FMX or DF.SNG only
- options to skip freemuxlet / use demuxlet
- decide which parameters go in `tool.config` vs `<sample>.config`
- possibly move the default QC cuts parameter to library-level, add use case where not present
- freemuxlet merge directory location?
