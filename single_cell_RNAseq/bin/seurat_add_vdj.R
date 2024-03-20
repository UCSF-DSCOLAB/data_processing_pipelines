library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)

args = commandArgs(trailingOnly=TRUE)
LIBRARY = args[1]
SOBJ = args[2]
DATA_TYPE = args[3] # TCR or BCR
CLONOTYPE_PATH = args[4]
CONTIG_PATH = args[5]
BASE_DIR = args[6]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

sobj = readRDS(SOBJ)

clonotype_data = load_clonotypes(LIBRARY, DATA_TYPE, CLONOTYPE_PATH, CONTIG_PATH)
sobj = AddMetaData(sobj, metadata=clonotype_data[colnames(sobj),])
saveRDS(sobj, file=sprintf("%s_w_%s.RDS", LIBRARY, DATA_TYPE), compress=F)