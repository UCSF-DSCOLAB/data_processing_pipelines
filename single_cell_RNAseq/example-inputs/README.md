
### Files
This `example-inputs` directory contains the following files:
* `default_qc_cuts.csv`: file with default QC cutoffs which can be helpful for visualizing initial QC. You will set the true cutoffs for each library, but this can provide a helpful starting point. Note that for particular tissues or datasets, you may want to adjust the default cutoffs for all libraries in a project to help better visualize between-library differences.
* `param_1.json`: an example parameter config file, which additionally shows how to specify a vcf and a variety of different data types. This parameter file does not refer to any existing toy data, it is just an example.
* `param_2.json`: example parameter config file for the toy dataset. See below for additional info on the toy dataset.
* `param_2v2.json`: updated version of `param_2.json` for `pre_qc` pipeline only.


### Toy Data
The toy data is available on c4 at:
`/krummellab/data1/pipeline_test_data/assays/scRNA_seq/modality/gex/downsampled_jurkat_tcell/inputs/fastqs/`

This full workflow test case was created from the a dataset of 50:50 jurkat-293t cells, originally from [Zheng et al 2017](https://www.nature.com/articles/ncomms14049), with full data available on the [10x website]( https://www.10xgenomics.com/resources/datasets/50-percent-50-percent-jurkat-293-t-cell-mixture-1-standard-1-1-0). A version of the bam that downsampled to 500 cells was downloaded from the [poscle tutorial](https://drive.google.com/drive/folders/1drNBY0SltMKpgLe_z9w1Swx1QK14uO5T), converted to fastqs. This was then artificially split into a dataset of multiple "pools" for testing purposes, with 200, 200, and 100 cells per library, and the first two libraries belonging to the same "pool". 

