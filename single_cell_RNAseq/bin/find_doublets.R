library(Seurat)
library(tidyverse)
library(DoubletFinder)

# Input arguments:
# 1. H5 file path
# 2. Freemuxlet samples file path (Use "NULL" if a single sample was loaded on 10x well)
# 3. Sample id/name (e.g. SICCA1-POOL-GN1-SCG1)
args = commandArgs(trailingOnly=TRUE)
CELLR_H5_PATH=args[1]
FMX_SAMPLE_PATH=args[2]

SAMPLE=args[3]
NCELLS_LOADED = as.numeric(args[4])
MINFEATURE=as.numeric(args[5])
MINCELL=as.numeric(args[6])

RANDOMSEED =args[7]
USE_INTER_DBL_RATE=args[8]
BASE_DIR=args[9]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

options(future.globals.maxSize= (2400*1024^2))


set.seed(RANDOMSEED)
NPCS_DOUBLETFINDER=30

print_message(sprintf(
  "
  -----------
  Running doublet finder with the following params:
    CELLR_H5_PATH=%s
    FMX_SAMPLE_PATH=%s
    SAMPLE=%s
    NCELLS_LOADED=%s
    MINFEATURE=%s
    MINCELL=%s
    RANDOMSEED=%s
    BASE_DIR=%s
  -----------
  ", CELLR_H5_PATH, FMX_SAMPLE_PATH, SAMPLE, NCELLS_LOADED, MINFEATURE, MINCELL, RANDOMSEED, BASE_DIR)
)


# adapted from Arjun's get10xMultipletRate.R
doublet_rate_df = data.frame(multiplet_rate_pct=c(0.40, 0.80, 1.60, 2.30, 3.10, 3.90, 4.60, 5.40, 6.10, 6.90, 7.60), 
                             num_cells_loaded=c(870, 1700, 3500, 5300, 7000, 8700, 10500, 12200, 14000, 15700, 17400),
                             num_cells_recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))

DBL_MODEL_recovered <- lm(multiplet_rate_pct ~ num_cells_recovered, data=doublet_rate_df)
DBL_MODEL_loaded <- lm(multiplet_rate_pct ~ num_cells_loaded, data=doublet_rate_df)
MODEL_recovered <- lm(num_cells_recovered ~ num_cells_loaded, data=doublet_rate_df)
MODEL_loaded <- lm(num_cells_loaded ~ num_cells_recovered, data=doublet_rate_df)

genDoubletTable = function(doublet_stats, ncells_loaded, nsamples){
  
  ncells = sum(doublet_stats$fmlDropletTypeComp)
  dbl_rate_10x_recovered = predict(DBL_MODEL_recovered, new=data.frame(num_cells_recovered=ncells)) / 100
  dbl_rate_10x_loaded = predict(DBL_MODEL_loaded, new=data.frame(num_cells_loaded=ncells_loaded)) / 100
  predicted_recovery = predict(MODEL_recovered, new=data.frame(num_cells_loaded=ncells_loaded))
  predicted_loaded = predict(MODEL_loaded, new=data.frame(num_cells_recovered=ncells))

  num_mx_dbls = ifelse("DBL" %in% rownames(doublet_stats$fmlDropletTypeComp), doublet_stats$fmlDropletTypeComp["DBL",][[1]], 10)
  num_mx_sngs = doublet_stats$fmlDropletTypeComp["SNG",][[1]]

  fmlDblRate = num_mx_dbls/ncells
  dblRateIntra = fmlDblRate/(nsamples-1)
  num_dbls_intra = dblRateIntra*ncells

  # doublet rate is fraction of *loaded* cells
  num_rem_dbl_recovered = dbl_rate_10x_recovered*predicted_loaded-num_mx_dbls 
  num_rem_dbl_loaded = dbl_rate_10x_loaded*ncells_loaded-num_mx_dbls 

  df_stat_tib = tibble(
    stat = c("ncells loaded", "ncells",  "predicted ncells recovered", "predicted loaded", "nsamples",
            "DBL rate (10x based on ncells recovered)", "DBL rate (10x based on ncells loaded)",
            "FMX interDBL rate", "FMX intraDBL rate", "Remaining doublets (10x recovered)", 
            "Remaining doublets (10x loaded)", "Remaining intra DBL (FMX)"),
    value = c(ncells_loaded, ncells, predicted_recovery, predicted_loaded, nsamples, 
      dbl_rate_10x_recovered, dbl_rate_10x_loaded, fmlDblRate, dblRateIntra, num_rem_dbl_recovered, 
      num_rem_dbl_loaded,  num_dbls_intra)
  )
  return(df_stat_tib)
}

runDoubletFinder <- function(sObj, freemuxlet=TRUE, use_inter_dbl_rate=TRUE, ncells_loaded) {
  effectiveIntraDblRate = 0 # Initialize
  nExp_poi=0
  if(freemuxlet & use_inter_dbl_rate ) {        # freemuxlet == true, means freemuxlet data is loaded to the meta.data slot.
    print_message("Using inter-individual doublet rate to estimate number of doublets")

    # remove filtered cells
    present.cells = rownames(sObj@meta.data[!is.na(sObj@meta.data$DROPLET.TYPE),])
    sngObj = subset(sObj, cells=present.cells)

    # freemuxlet doublet (DBL) rate
    fmlDblRate = ifelse("DBL" %in% rownames(sngObj@misc$scStat$fmlDropletTypeProp), 
      sngObj@misc$scStat$fmlDropletTypeProp["DBL",], 10)

    # Determine the number of samples pooled per 10x lane from the freemuxlet output.
    
    ## TODO: do this based on config?
    noOfSmpls = length( unique( sngObj$BEST.GUESS[ sngObj$DROPLET.TYPE == "SNG" ] ) )

    dblRateIntra = fmlDblRate/(noOfSmpls-1)

    # Calculate adjustment factor:
    #  dblRateIntra gives the rate of intra-ind doublets out of the total
    #  but we will remove inter-ind DBL *BEFORE* running doubletFinder
    #  so we need to adjust for the change in denominator
    num_singlets = sngObj@misc$scStat$fmlDropletTypeComp[c("SNG"),]
    #num_singlets_and_doublets = sum(sngObj@misc$scStat$fmlDropletTypeComp[c("SNG","DBL"),] )
    adjustFactor = ncol(sngObj)/num_singlets

    # Effective doublet rate after accounting for removal of DBL cells.
    effectiveDblRate = dblRateIntra*adjustFactor 

    sngObj = subset(sngObj, DROPLET.TYPE=="SNG")
    nExp_poi <- ceiling( effectiveDblRate * ncol(sngObj))# remaining doublets
  } else {
    if (freemuxlet & !use_inter_dbl_rate){
      print_message("Using cells recovered to estimate number of doublets")
      present.cells = rownames(sObj@meta.data[!is.na(sObj@meta.data$DROPLET.TYPE),])
      sngObj = subset(sObj, cells=present.cells)
      num_inter_dbl = sum(sngObj$DROPLET.TYPE=="DBL")
      sngObj = subset(sngObj, DROPLET.TYPE=="SNG")
      effectiveDblRate = unname(predict(DBL_MODEL_recovered, new=data.frame(num_cells_recovered=ncol(sngObj))) / 100)
      predicted_loaded = unname(predict(MODEL_loaded, new=data.frame(num_cells_recovered=ncol(sngObj))))
      nExp_poi = ceiling(predicted_loaded*effectiveDblRate) - num_inter_dbl
     } else {
       print_message("No freemuxlet data, using number of cells loaded to estimate")
       sngObj = sObj
       effectiveDblRate = unname(predict(DBL_MODEL_loaded, new=data.frame(num_cells_loaded=ncells_loaded)) / 100)
       nExp_poi = ceiling(ncells_loaded*effectiveDblRate) 
     } 

  }

  names(effectiveDblRate) = NULL
  sObj@misc$scStat$effDblRate = effectiveDblRate
  print_message("Effective intra-sample doublet rate: ", round(effectiveDblRate * 100, 4), "%")
  print_message("Number of remaining doublets: ", nExp_poi)

  # Now that the effective intra-sample doublet rate has been calculated, let's run doubletFinder 
  doubletF_stat = list()
  doubletFinderCalls = c()
 
  sng.cells=colnames(sngObj)
  doubletF_stat[["nCells"]] = length(sng.cells)

  # If the number of cells is smaller than npcs, adjust the npcs.
  npcs = NPCS_DOUBLETFINDER  # Setting number of PCs/dimensions to use for UMAP/Clustering.
  if(length(sng.cells) <= npcs) npcs = length(sng.cells) - 1
  print_message("\tThe number of PCs used for DoubletFinder analysis: ", npcs)
  doubletF_stat[["nPCs"]] = npcs

  print_message("\tRunning Seurat pipeline to make SCT-based UMAP/clustering for the use in DoubletFinder")
  sngObj <- SCTransform(sngObj, verbose = FALSE)
  sngObj <- RunPCA(sngObj, npcs=npcs, verbose = FALSE)
  sngObj <- RunUMAP(sngObj, dims = 1:npcs, verbose = FALSE)
  sngObj <- FindNeighbors(sngObj, dims = 1:npcs, verbose = FALSE)
  sngObj <- FindClusters(sngObj, verbose = FALSE)
  
  print_message("\tEstimating pK")
  sweep.res.list <- paramSweep_v3(sngObj, PCs = 1:npcs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK = as.numeric(as.vector(bcmvn[ bcmvn$BCmetric == max(bcmvn$BCmetric) ,]$pK))
  print_message("\tChosen pK: ", pK)
  doubletF_stat[["pK"]] = pK

  # The modeling of homotypic doublets requires celltype annotation. In place of the celltype annotation, I am using the cluster identity as celltype labels (as suggested by the authors of doubletFinder) - since clusters should correspond to different celltypes anyway.
  annotations = sngObj@active.ident
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sObj@meta.data$ClusteringResults

  nExp_poi.adj <- ceiling(nExp_poi*(1-homotypic.prop))
  print_message("\tEstimated no. of heterotypic doublets:", nExp_poi.adj)

  print_message("\tRunning DoubletFinder")
  sngObj <- doubletFinder_v3(sngObj, PCs = 1:npcs, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  score=grep("^pANN", colnames(sngObj@meta.data))[1]
  ind = grep("^DF", colnames(sngObj@meta.data))[1]
  dfCalls = data.frame(sngObj@meta.data[,c(ind, score)])
  colnames(dfCalls)=c("doubletFinderCalls", "doubletFinderScore")
  rownames(dfCalls) = rownames(sngObj@meta.data)
  doubletF_stat[["nDFSng"]] = sum(dfCalls$doubletFinderCalls == "Singlet")
  doubletF_stat[["nDFDbl"]] = sum(dfCalls$doubletFinderCalls == "Doublet")
  sObj@misc$scStat$doubletF_stat = doubletF_stat

  sObj = AddMetaData(sObj, metadata=dfCalls)

  return(sObj)
}


data = Read10X_h5(CELLR_H5_PATH)
if ("Gene Expression" %in% names(data)){
  seuratObj = CreateSeuratObject(data$`Gene Expression`, project = SAMPLE, 
                                 min.cells = MINCELL, 
                                 min.features = MINFEATURE)
} else {
  seuratObj = CreateSeuratObject(data, project = SAMPLE, 
                                 min.cells = MINCELL, 
                                 min.features = MINFEATURE)
}
print_message("Dimension of the GEX count data:", paste(dim(seuratObj), collapse=" x ") )
seuratObj@misc$scStat = list()

# Load freemuxlet calls
use_inter_dbl_rate = (USE_INTER_DBL_RATE=="true")
freemuxlet = TRUE
if(is.null(FMX_SAMPLE_PATH)) {
  freemuxlet = FALSE
  use_inter_dbl_rate = FALSE
} else {
  seuratObj = loadFreemuxletData(seuratObj, FMX_SAMPLE_PATH)
  # TODO: should this cause an error?
  if (!"SNG" %in% rownames(seuratObj@misc$scStat$fmlDropletTypeComp)){
    print_message("Warning - there are no singlets in your demultiplexing output. 
                  Please re-run free/demuxlet with alternate parameters or filtering before running doublet finder.")
    exit()
  }
}





# Identify intra-sample doublets
seuratObj = runDoubletFinder(seuratObj, freemuxlet, use_inter_dbl_rate, NCELLS_LOADED)
seuratObj@meta.data$DROPLET.TYPE.FINAL = ifelse(
    seuratObj$DROPLET.TYPE == "AMB", "AMB", ifelse(
      seuratObj$DROPLET.TYPE == "DBL", "Inter.DBL", ifelse(
        seuratObj$doubletFinderCalls == "Doublet", "Intra.DBL", "SNG"
      ) 
    )
  )

tbl = table(seuratObj$DROPLET.TYPE.FINAL)
seuratObj@misc$scStat$finalDropletTypeComp = as.matrix(tbl)
seuratObj@misc$scStat$finalDropletTypeProp = as.matrix(prop.table(tbl))

saveRDS(seuratObj@misc$scStat, file=sprintf("%s_doubletStat.rds", SAMPLE))
saveRDS(seuratObj, file=sprintf("%s_seurat_object_findingDoublets.rds", SAMPLE))

nsamples = length( unique( seuratObj$BEST.GUESS[ seuratObj$DROPLET.TYPE == "SNG" & !is.na(seuratObj$DROPLET.TYPE) ] ) )

doublet_tab = genDoubletTable(seuratObj@misc$scStat, NCELLS_LOADED, nsamples)
doublet_tab %>%
  mutate(value=round(value, 3)) %>% 
  write_csv(sprintf("%s_doublet_counts.csv", SAMPLE))