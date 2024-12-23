
This code was initially developed at https://github.com/UCSF-DSCOLAB/sc_seq_nextflow/. Please note that the original
repository is no longer maintained and has been marked as deprecated. Going forward, all enhancements, bug fixes, 
and active development are to be carried out in this repository.

## Summary

This pipeline is a nextflow implementation of the `sicca_sc-seq_pipeline` written by Kim Taylor and Ravi Patel, with the following modifications: (1) addition of a freemuxlet merge step to combine across pools and (2) addition of an option to filter based on QC pre-freemuxlet. Future efforts will make this more flexible for different inputs.

For non-tissue data, Steps 1-3 are performed with `pre_qc` and step 5 is performed with `post_qc`

The pipeline performs the following steps:
1. Alignment with cellranger (optional - you can start with aligned data)
2. Freemuxlet or Demuxlet: Deconvolute sample identities and finds interindividual doublets (unsupervised vs. supervised) 
3. DoubletFinder: Used to identify intraindividual doublets (optional). This is determined by the number of inter-individual doublets from freemuxlet/demuxlet.
4. Filtering: Manually set qc filters (Note: this can be performed before step 2 with `pre_fmx_qc` and `post_fmx_qc` step path, which is recommended for tissue data with high ambient RNA. Further details below.)
5. Seurat: Normalize and scale RNA data, perform dimensionality reduction, and generate initial clusterings and visualizations of the data. (ADT data will be DSB-normalized)


## How to run on the c4 cluster

Although this pipeline is not specific to the HPC cluster we use within DSCoLab, we have not tested it outside of c4.

### Pre-reqs: Java and Nextflow
You can add the following to your `~/.bashrc`, then run `source ~/.bashrc` to connect with publicly-accessible locations where we have installed nextflow and java:

```bash
# global environment variables
export JAVA_HOME="/krummellab/data1/software/jdk/jdk-17.0.5"
export NXF_JAVA_HOME="/krummellab/data1/software/jdk/jdk-17.0.5"
export PATH=$PATH:"/krummellab/data1/software/nextflow/22.10.4_build_5836/"
```

### Before running
First, you will need to create a config file that is adjusted for your target data.
For example config files, see [example-inputs/param_1.json](example-inputs/param_1.json) and [example-inputs/param_2.json](example-inputs/param_2.json). Additional documentation is provided in [example-inputs/README.md](example-inputs/README.md).


Within these files, you will need to update certain parameters in order to specify:
1. "project_dir" should point to the /krummellab/data1/immunox/<project_name> directory for your project. Standard DSCoLab structure for the contents of this directory are expected: Fastqs are expected to come from within 'data/<modality>/raw' folders within here, and outputs will be created within 'data/single_cell_GEX/processed/<pool_names>' here as well.
2. any "settings" should be adjusted for your data type and how you want your run to be configured.
3. "pools" contents must be adjusted to give the names and characteristics of your own libraries

### How to run

**To run:**
`sbatch run.sh <path_to_config.json> <step>`

&nbsp;&nbsp;&nbsp;&nbsp; A more literal example: `sbatch run.sh path/to/my/config.json pre_qc`

&nbsp;&nbsp;&nbsp;&nbsp; Also be sure to take note of the number reported back to you as "Submitted batch job #######" after running this command. That number is what to use as <job_id_for_original_run> if you have need to resume the job later.

**To resume a previous run:**
 `sbatch run_resume.sh <path_to_config.json> <step> <job_id_for_original_run>`

**The `step` input:**

`step` must be one of ['pre_qc', 'post_qc', pre_fmx_qc', 'post_fmx_qc'].
The pipeline is designed to be run in either of 2 ways, and both involve an internal breakpoint where you will need to supply cell filtration cutoffs:
- Path 1 is the "standard" method, run with (1) step as `pre_qc`, (2) setting of cutoffs, then (3) step as `post_qc`.
- Path 2 allows for QC filtering to be performed before freemuxlet/demuxlet which can be useful for tissue data, and is used with: (1) step as `pre_fmx_qc`, (2) setting of cutoffs, then (3) step as `post_fmx_qc`.

(Note: we have also added a `df_auto` step to run doubletfinder and automating processing only. This is designed to be run only in special circumstances where the user decides to re-run only these steps. If you have run the pre_fmx path and reviewed cutoffs, it will complete. Otherwise, it will fail at the cutoff reviewing stage and you will need to complete with `post_qc`.)

### After run cleanup
By default, the nextflow working directory will be:
`/c4/scratch/<user>/nextflow/<original_job_id>`
This directory is automatically deleted on successful completion of an initial run, thus there is often no manual cleanup needed.

However, if the pipeline fails, the directory is left in place. Please deleted this directory if you do not plan to try again, or you can resume a failed run with the `run_resume.sh` method described in the "How to run" section. When resuming a run this way, the same working directory will be used, but regardless of whether the resumed run is successful or not the directory will not be deleted.
**After your resumed run completes, or you decide not to follow up a failed run, please manually remove this directory.**

### Example run with toy data
0. Download and unzip this github.
1. Create a directory for this run, let's call it `${TOY_PROJECT_DIR}`
2. Create the following subdirectiories:
`mkdir -p ${TOY_PROJECT_DIR}/data/single_cell_GEX/raw/; mkdir -p ${TOY_PROJECT_DIR}/freemuxlet_data/`
3. Copy the toy data from c4 to your raw directory `cp /krummellab/data1/pipeline_test_data/assays/scRNA_seq/modality/gex/downsampled_jurkat_tcell/inputs/ ${TOY_PROJECT_DIR}/data/single_cell_GEX/raw/`
4. Make a copy of the run's config file.
`cp example-inputs/param_2.json example-inputs/my_toy_config.json`

5. Open the run's config file and edit the first two directiories to point to `${TOY_PROJECT_DIR}`, and `${TOY_PROJECT_DIR}/freemuxlet_data/` respectively. The third should point to the directory this repository is in, specifically the subdirectory: `${DATA_PROCESSING_PIPELINE_REPO}/single_cell_RNAseq/example-inputs/`
6. Submit the run. Note we are using `-profile test` for these test data because they are much smaller. Be sure to remove this flag for any real run as it scales down the resources requested for each task to levels that are not viable for real data. Note that a test run will use the freecycle partition (see more on partitions in additional details).
`sbatch run.sh example-inputs/my_toy_config.json pre_qc -profile test`
7. Once complete, manually examine the diagnostic plots and set cutoffs. Then, run the final step: `sbatch run.sh example-inputs/my_toy_config.json post_qc -profile test`

## Additional Details

### Notes and limitations

* the pipeline assumes the standard immunox directory structure
* data *must* be pooled libraries that require genetic demultiplexing to use the pipeline
* the references are all for human hg38, you will need to edit `config/reference.config` if you want different refs or have non-human data
* this branch does not yet work for VDJ (this functionality will be added soon)
* this does not work for HTO data
* for snRNA-seq: you will need to edit the cellranger step to include introns if you are using cellranger < v7.0.0

### Side-effects of the Default QC Cutoffs
In addition to being used as the initial values for your 'qc_cuts.csv', the values in the file pointed to by your config file's "default_qc_cuts_dir" & "default_qc_cuts_file" elements also determine where initial lines are drawn in QC plots output by the `pre_qc` and `pre_fmx_qc` steps.

Of note, it can be useful to adjust these values sensibly before / for all runs by adjusting this file.  This can be especially useful with tissue data where it is often useful to compare across libraries and batches, with consistent cutoff values, as a method of assessing relative quality across the entire project.

### c4 partitions, additional configuration
We have the following partitions available on the c4 cluster for DSCoLab members: krummellab,common (default partitions, which are set in the environment variable `$SBATCH_PARTITION`) and freecycle. Freecycle allows us to run jobs on any compute that is not in use (24h job max), with the caveat that the job may be cancelled if a partition owner requires these resources. We have set the default for `test` jobs to freecycle, and for standard jobs to `freecycle,krummellab,common`, which means they are first submitted to freecycle, then krummellab, then common. If you do not have access to the krummellab partition on c4, you may to change the process.queue parameter in your nextflow.config file to reflect the partitions you have access to, otherwise it will default to freecycle and common. You may want to remove freecycle and/or other cluster options, such as the errorStrategy (what happens to the whole pipeline when it errors). See comments in `nextflow.config` describing what we have set as defaults, and what you may want to adjust.



### Data Generation

For those that are not familiar with sequencing and data generation, recall that:

1. We obtain $s$ biological samples / biospecimens. Each sample maps to a unique individual and  
contains $c$ cells.
2. For each sample we isolate $n$ cells, and then load those cells into the sequencer. When $s > 1$ this is
called pooling. 
3. The sequencer runs and generates a library. Each library can contain $r$ reads.

Note: Pooling is somewhat specific to Colabs

### What else might I need to configure?

There are multiple less standard inputs that are required for pipelines to run, including:

- Location of reference genomes and containers
- Process specific settings and flags

In order to view what all of these settings are, you can check out `nextflow.config`.
To actually supply the parameters to this pipeline, you must submit a json file with these values.
Some examples include: `example-inputs/param_1.json` `example-inputs/param_2.json`.

There is also a directory called `config/` however these are settings that specific to c4 and typically
do not need to be tweaked. `nextflow.config` imports these settings for you.

---

## Further notes specific to future developers
### Testing

#### Regression tests

You will need nf-test installed. Please follow the instructions here: https://github.com/askimed/nf-test#installation

To run regression tests: 

`cd single_cell_RNAseq`

`/krummellab/data1/software/bin/nf-testtest tests/pipeline_pre_qc.nf.test --without-trace`

#### How do the tests work

- For each test we copy all the necessary file inputs from `/krummellab/data1/integration_test_user/tutorial_lib_sep` 
into an isolated testing environment.
- That testing environment is `.nf-test/hash/meta`
- We inject these input paths into our pipeline.
- The tests then operate on these files.
- We assert that the necessary files are created in this testing environment.

#### Test data

In general test data should be located in: `/krummellab/data1/pipeline_test_data/`.

Currently we only have GEX (Gene expression data) to use as a means of testing our pipeline.

The fastqs and variants can be found in the following directory:

`/krummellab/data1/pipeline_test_data/assays/scRNA_seq/modality/gex/downsampled_jurkat_tcell/inputs/`

We also want:

- CITE (Protein level data)
- BCR (B-cell receptor)
- TCR (T-cell receptor)


### Conventions


#### Params

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

#### Variable names

Groovy is JVM-esque so it might have made sense to use `camelCase` for variable names.
However, most of the code written in this repo is `snake_case`, so lets continue with this
convention.

It is also helpful to prefix channels with the `ch_` prefix.



