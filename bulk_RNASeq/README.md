# RNAX
A computational pipeline for variant calling and expression analysis of bulk RNA-sequencing data.

N.B: The pipeline is under active development and not yet suitable for deployment.

## Installation (for development only)
* Create virtual environment using Anaconda
```bash
conda env create -f environment.yml -n rnax_env
```

## Usage (for development only)
* Create sample sheet specifying each sample, locations of their sequencing reads (FASTQ), and whether they are single-end or paired-end reads. The figure below is an example of a correctly formatted sample sheet:
![sample_sheet](docs/figs/sample_sheet_example.png)
* Run the DSL2 pipeline
```bash
nextflow run bulk_rna_seq.nf -w [working_directory]
```

## Authors
Emily Flynn
Al Latif
Daniel Bunis
Walter Eckalbar