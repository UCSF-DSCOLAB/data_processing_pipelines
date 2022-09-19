
7/28/2022
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz


wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt




Checked the md5sums


PARENT_DIR=/krummellab/data1/erflynn/ref/analysis_set/
SOFTWARE_DIR=/krummellab/data1/erflynn/software/
# containers used + versions
SINGULARITY_DIR=/krummellab/data1/singularity_images/
STAR_CONTAINER=${SINGULARITY_DIR}/STAR/2.7.8a/STAR.sif
RSEM_CONTAINER=${SINGULARITY_DIR}/RSEM/1.3.3/RSEM.sif
KALLISTO_CONTAINER=${SINGULARITY_DIR}/kallisto/0.46.1/kallisto.sif
SAMTOOLS_CONTAINER=${SINGULARITY_DIR}/samtools/1.12/samtools.sif
PICARD_CONTAINER=${SINGULARITY_DIR}/picard/2.25.3/picard.sif

cd ${PARENT_DIR}
mkdir star
mkdir rsem
mkdir kallisto


gunzip hg38.analysisSet.fa.gz 

gunzip hg38.knownGene.gtf.gz 

# 1. set up STAR index

singularity exec -B ${PARENT_DIR} \
    ${STAR_CONTAINER} \
  STAR --runThreadN 11 \
 --runMode genomeGenerate \
 --genomeDir ${PARENT_DIR}/star_chr/ \
 --genomeFastaFiles ${PARENT_DIR}/hg38.analysisSet.fa \
 --sjdbGTFfile ${PARENT_DIR}/hg38.knownGene.gtf
 --sjdbOverhang 149 # readLength -1


# 2. set up RSEM index
singularity exec -B ${PARENT_DIR} \
    ${RSEM_CONTAINER} \
    rsem-prepare-reference --gtf ${PARENT_DIR}/hg38.knownGene.gtf \
               ${PARENT_DIR}/hg38.analysisSet.fa \
               ${PARENT_DIR}/rsem/hg38

# 3. Create a flat file
${SOFTWARE_DIR}/gtfToGenePred -genePredExt ${PARENT_DIR}/hg38.knownGene.gtf ${PARENT_DIR}/hg38.knownGene_flat0.txt

awk -v OFS='\t' '{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' hg38.knownGene_flat0.txt > hg38.knownGene_flat.txt


# 3b. exon only files

# TODO: filter to remove contigs that are on not in the dict/transcriptome

cat hg38.knownGene.gtf | awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > exon_only0.bed

grep -v "_alt" exon_only0.bed | grep -v "_fix" >  exon_only.bed

cat hg38.knownGene.gtf | awk '{if ($3=="exon") {print $0}}'  > exon_only.gtf


# 4. create ribosome file
singularity exec -B $PARENT_DIR \
            ${SAMTOOLS_CONTAINER} \
            samtools faidx ${PARENT_DIR}/hg38.analysisSet.fa

cat hg38.analysisSet.fa.fai | awk 'BEGIN {print "@HD\tVN:1.0\tSO:coordinate"} {printf "@SQ\tSN:%s\tLN:%s\n", $1, $2}' > ribo_intervals.txt

grep "transcript_type \"rRNA" hg38.knownGene.gtf | awk -v OFS='\t' '{print $1, $4, $5, $7, $16}' | uniq | sed -E 's/(\"|;)//g'  >> ribo_intervals.txt 


# 5. Create kallisto index - first have to generate tx file
${SOFTWARE_DIR}/gffread/gffread-0.12.7/gffread -w ${PARENT_DIR}/hg38.transcriptome.fa -g ${PARENT_DIR}/hg38.analysisSet.fa ${PARENT_DIR}/hg38.knownGene.gtf 


#NOTE: RESULTED IN THIS ERROR:
# Warning: couldn't find fasta record for 'chr10_GL383545v1_alt'!
#Error: no genomic sequence available (check -g option!).

singularity exec -B ${PARENT_DIR} ${KALLISTO_CONTAINER} \
            kallisto index -i ${PARENT_DIR}/kallisto/hg38.index \
             --make-unique ${PARENT_DIR}/hg38.transcriptome.fa

# 6. Sequence dictionary
singularity exec -B ${PARENT_DIR} ${PICARD_CONTAINER} java -jar /opt/picard/picard.jar CreateSequenceDictionary R=${PARENT_DIR}/hg38.analysisSet.fa  O=${PARENT_DIR}/hg38.analysisSet.dict


Contig chr10_GL383545v1_alt given as location, but this contig isn't present in the Fasta sequence dictionary
exon_only.bed

GCF_000001405.39_mapped.txt.gz

