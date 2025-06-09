# Reference Genome Preparation for Bulk RNAseq processing
A computational pipeline for building and preparing reference genome files for the associated computational pipeline for Bulk RNAseq processing.

## Installation
* Clone repository
```bash
git clone git@github.com:UCSF-DSCOLAB/data_processing_pipelines.git
```
* Install a modern JAVA version
```bash
wget https://download.oracle.com/java/17/archive/jdk-17.0.7_linux-x64_bin.tar.gz
```
* Decompress the downloaded file
```bash
tar -xzvf jdk-17.0.7_linux-x64_bin.tar.gz
```
* Update your `~/.bashrc` file to point to the modern JDK
```bash
# Add global environment variables for JDK
export JAVA_HOME="/path/to/jdk-17.0.6"
export NXF_JAVA_HOME="path/to/jdk-17.0.6"
```
* Download and install Nextflow
```bash
wget -qO- https://get.nextflow.io | bash
```
* Make Nextflow an executable
```bash
chmod +x nextflow
```

## Usage
* Load the Anaconda environment manager
```bash
module load CBI miniconda3/23.3.1-0-py39
```
* Open the `config/nextflow.config` file and specify your file paths
    * `genome`: path to the reference genome FASTA file
    * `gtf`: path to the reference genome GTF file
    * `tmp_dir`: path to the folder for storing temporary files e.g. /scratch/username
    * `index_vcf`: an optional parameter for generating an index to the reference VCF file
    * `dbsnp`: path to the reference genome VCF file (optional, can be empty with `index_vcf` set to `false`)
    * `reference_directory`: path to the folder for outputting the reference files
* N.B: Remember to save the changes to the file after editing
* Run the DSL2 pipeline to build the reference genome
```
run_prep_c4.sh
```