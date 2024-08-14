
### Files
This `example-inputs` directory contains the following files:
* `default_qc_cuts.csv`: file with default QC cutoffs which can be helpful for visualizing initial QC. You will set the true cutoffs for each library, but this can provide a helpful starting point. Note that for particular tissues or datasets, you may want to adjust the default cutoffs for all libraries in a project to help better visualize between-library differences.
* `param_1.json`: an example parameter config file for the `pre` and `post` pipelines. Which additionally shows how to specify a vcf and a variety of different data types. This parameter file does not refer to any existing toy data, it is just an example. 
* `param_2.json`: an example parameter config file for the toy dataset (`pre` and `post` pipelines). See below for additional info on the toy dataset.
* `fmx_param_1.json`: Similar to `param_1.json` except its format is intended for the `pre_fmx` and `post_fmx` pipelines.


### Toy Data
The toy data is available on c4 at:
`/krummellab/data1/pipeline_test_data/assays/scRNA_seq/modality/gex/downsampled_jurkat_tcell/inputs/fastqs/`
`/krummellab/data1/pipeline_test_data/assays/scRNA_seq/modality/gex/downsampled_jurkat_tcell/inputs/variants/`

This full workflow test case was created from the a dataset of 50:50 jurkat-293t cells, originally from [Zheng et al 2017](https://www.nature.com/articles/ncomms14049), with full data available on the [10x website]( https://www.10xgenomics.com/resources/datasets/50-percent-50-percent-jurkat-293-t-cell-mixture-1-standard-1-1-0). A version of the bam that downsampled to 500 cells was downloaded from the [poscle tutorial](https://drive.google.com/drive/folders/1drNBY0SltMKpgLe_z9w1Swx1QK14uO5T), converted to fastqs. This was then artificially split into a dataset of multiple "pools" for testing purposes, with 200, 200, and 100 cells per library, and the first two libraries belonging to the same "pool". The corresponding VCF was downloaded from Box [https://ucsf.app.box.com/s/vg1bycvsjgyg63gkqsputprq5rxzjl6k], a chromosome prefix was then added, the results lifted over from hg19 to hg38 using Picard tools.

### Parameter file fields explained

Below is the json for `param_2.json`, with comments added to describe what each line does.
Note that comments cannot be included in json used for a run.

```
{ "project_dir" : "tutorial_lip_sep", # the location of the project directory. your raw files must be in `data/single_cell_<data_type>/raw/<library>`, and all outputs will be written to this directory in `data/single_cell_<data_type>/processed/<library>`
    "settings" : {
    "add_tcr" : false, # ignore this, functionality not yet implemented
    "add_bcr" : false, # ignore this, functionality not yet implemented
    "skip_cellranger": false, # whether or not to skip cellranger -- set this to true if you have already aligned the data to avoid re-running a computationally intensive step. if you set this to true, the processed cellranger output must be located in `data/single_cell_<data_type>/processed/<library>/cellranger/` for the subsequent steps to work
    "merge_for_demux" : true, # whether to merge across libraries in pool prior to running free/demuxlet. this should improve the demultiplexing by providing additional data to the algorithm
    "merge_demux_dir" : "freemuxlet_data/", # where to put the merged data for free/demuxlet. these data will be unmerged and stored in the appropriate `processed/<library>` directory, but it can be helpful to keep the merged data for future re-runs
    "demux_method" : "freemuxlet", # this can be one of freemuxlet or demuxlet
    "fmx_assign_to_gt": false, # whether to assign to genotypes
    "ref_vcf_type": "bulk", # can be "bulk" or "array", the source of the vcf for gtcheck
    "run_doubletfinder" : true, # whether to run doubletfinder to ID intra-individual doublets. the number of intra-individual doublets will be based on the number of inter-individual doublets identified by free/demuxlet
    "mincell" : 3, # parameter for load into seurat -- minimum number of cells required to keep a gene
    "minfeature" : 100, # parameter for load into seurats - minimum number of features required to keep a cell
    "default_qc_cuts_dir": "example-inputs", # location of default qc cuts file -- this can be the path to repo 
    "default_qc_cuts_file": "default_qc_cuts.csv", # default qc cuts file, see notes above about this
    "randomseed" : 21212, # random seed for reproducibility
    "remove_demux_DBL": true, # whether to remove free/demux doublets
    "remove_all_DBL": true	# whether to remove all doublets (free/demux + doubletfinder)
  },
  "pools" : [
    { "name" : "DM1", # name of the pool
      "nsamples" : "2", # number of samples in the pool, needed for freemuxlet
      "vcf": "/path_to/jurkat_293t_exons_only_w_chr_hg38.vcf ", # vcf containing just the individuals in the pool, needed for demuxlet or fmx_assign_to_gt, can be left blank otherwise. this must be the entire path
      "libraries": [
        { "name": "TEST-POOL-DM1-SCG1",  # name of the library
          "ncells_loaded": 200, # number of cells loaded, helpful for checking counts and relative doublets
          "data_types": ["GEX"]  # should be a list [], must contain one of GEX or CITE (cannot contain both), and then additional modalities such as TCR or BCR
        },
        {"name" : "TEST-POOL-DM1-SCG2",
          "ncells_loaded": 200,
          "data_types": ["GEX", "TCR"]
        }
      ]
    },
    { "name": "DM2" ,
        "nsamples" : "2",
        "vcf": "",
        "libraries": [
          {"name": "TEST-POOL-DM2-SCG1",
            "ncells_loaded": 100,
            "data_types": ["CITE", "BCR", "TCR"]
          }
        ]
    }
  ]
}
```
