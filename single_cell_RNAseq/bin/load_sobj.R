library(Seurat)
library(tidyverse)

# Input arguments:
# 1. H5 file path
# 2. Freemuxlet samples file path (Use "NULL" if a single sample was loaded on 10x well)
# 3. Sample id/name (e.g. SICCA1-POOL-GN1-SCG1)
args = commandArgs(trailingOnly=TRUE)
CELLR_H5_PATH=args[1]
FMX_SAMPLE_PATH=args[2]
SAMPLE=args[3]
MINFEATURE=as.numeric(args[4])
MINCELL=as.numeric(args[5])
BASE_DIR=args[6]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

print_message(sprintf(
  "
  -----------
  Running seurat with the following params:
    CELLR_H5_PATH=%s
    FMX_SAMPLE_PATH=%s
    SAMPLE=%s
    MINFEATURE=%s
    MINCELL=%s
    BASE_DIR=%s
  -----------
  ", CELLR_H5_PATH, FMX_SAMPLE_PATH, SAMPLE, MINFEATURE, MINCELL, BASE_DIR)
)



data = Read10X_h5(CELLR_H5_PATH)
if ("Gene Expression" %in% names(data)){
  seuratObj = CreateSeuratObject(data$`Gene Expression`, project = SAMPLE, 
                                 min.cells = MINCELL, 
                                 min.features = MINFEATURE)
  print(dim(data$`Gene Expression`))
} else {
  seuratObj = CreateSeuratObject(data, project = SAMPLE, 
                                 min.cells = MINCELL, 
                                 min.features = MINFEATURE)
}


print_message("Dimension of the GEX count data:", paste(dim(seuratObj), collapse=" x ") )
seuratObj@misc$scStat = list()

# Load freemuxlet calls
freemuxlet = TRUE
if(FMX_SAMPLE_PATH=="null") {
  freemuxlet = FALSE
} else {
  seuratObj = loadFreemuxletData(seuratObj, FMX_SAMPLE_PATH)
}
saveRDS_(seuratObj, file=sprintf("%s_initial_raw.rds", SAMPLE))