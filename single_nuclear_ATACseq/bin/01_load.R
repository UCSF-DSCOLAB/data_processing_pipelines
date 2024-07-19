library(ArchR)
library(tidyverse)
library(yaml)
set.seed(1)

PROJECT="my_project" # fill in your project

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR=sprintf("/krummellab/data1/immunox/%s/data/single_nuclear_ATAC/processed", PROJECT) 
WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)
OUT_DIR=sprintf("%s/%s/cell_filter/", PARENT_DIR, LIBRARY)
dir.create(WORKING_DIR, showWarnings=T)
dir.create(OUT_DIR, showWarnings=T)

addArchRThreads(threads = 4) 
addArchRGenome("hg38")

DEFAULT_TSS_MIN=4
DEFAULT_NFRAG_MIN=1000

list_inputs = ifelse(file.exists(sprintf('%s/%s/cellranger/atac_fragments.tsv.gz', PARENT_DIR, LIBRARY)),
    sprintf('%s/%s/cellranger/atac_fragments.tsv.gz', PARENT_DIR, LIBRARY),
    sprintf('%s/%s/cellranger/fragments.tsv.gz', PARENT_DIR, LIBRARY))

setwd(WORKING_DIR)

ArrowFiles <- createArrowFiles(
  inputFiles = list_inputs,
  sampleNames = LIBRARY,
  minTSS = DEFAULT_TSS_MIN,
  minFrags = DEFAULT_NFRAG_MIN, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# write out the cutoffs to a YAML file
cutoffs = list("tss.min"= DEFAULT_TSS_MIN, "nFrag.min"=DEFAULT_NFRAG_MIN, "reviewed"="FALSE")

write_yaml(cutoffs, sprintf("%s/cutoffs.yml", OUT_DIR))

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sprintf("%s/results", WORKING_DIR),
  copyArrows = TRUE 
)
saveArchRProject(proj)


df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

df %>%
  as_tibble(rownames="BARCODE") %>%
  separate(BARCODE, into=c("LIBRARY", "cell_id"), sep="#", remove=F) %>%
  write_csv(file=sprintf("%s/cell_qc_meta.csv", OUT_DIR))

p= ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = T,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))) + 
  geom_hline(yintercept = DEFAULT_TSS_MIN, lty = "dashed") + 
  geom_vline(xintercept = log10(DEFAULT_NFRAG_MIN), lty = "dashed")
pdf(sprintf("%s/unfilt_tss_vs_nFrag.pdf", OUT_DIR), height=5, width=5)
p
dev.off()


