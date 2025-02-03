library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)

args = commandArgs(trailingOnly=TRUE)
LIBRARY = args[1]
SOBJ = args[2]
DATA_TYPE = args[3] # TCR or BCR
CONTIG_PATH = args[4]
BASE_DIR = args[5]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

print_message(sprintf(
  "
  -----------
  Running seurat add vdj with the following params:
    LIBRARY=%s
    SOBJ=%s
    DATA_TYPE=%s
    BASE_DIR=%s
    CONTIG_PATH=%s
  -----------
  ", LIBRARY, SOBJ, DATA_TYPE, BASE_DIR, CONTIG_PATH)
)


sobj = readRDS(SOBJ)

clonotype_data = load_clonotypes(LIBRARY, DATA_TYPE, CONTIG_PATH)
sobj = AddMetaData(sobj, metadata=clonotype_data[colnames(sobj),])
saveRDS(sobj, file=sprintf("%s_w_%s.RDS", LIBRARY, DATA_TYPE), compress=F)