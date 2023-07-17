library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)

args = commandArgs(trailingOnly=TRUE)
LIBRARY = args[1]
SOBJ = args[2]
DATA_TYPE = args[3] # TCR or BCR
BASE_DIR=args[4]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

sobj = readRDS(SOBJ)

clonotype_data = load_clonotypes(LIBRARY, DATA_TYPE)
sobj = AddMetaData(sobj, metadata=clonotype_data)
saveRDS(sobj, file=sprintf("sobj_w_%s.RDS", DATA_TYPE), compress=F)