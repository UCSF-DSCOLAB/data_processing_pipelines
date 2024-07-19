library(ArchR)
library(tidyverse)
set.seed(1)

PROJECT="my_project" # FILL IN
args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR=sprintf("/krummellab/data1/immunox/%s/data/single_nuclear_ATAC/processed", PROJECT) 

WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)
dir.create(WORKING_DIR)
#AMULET_DIR=sprintf("%s/%s/amulet/", PARENT_DIR, LIBRARY)
FILTER_DIR=sprintf("%s/%s/cell_filter/", PARENT_DIR, LIBRARY)

setwd(WORKING_DIR)

# load
proj = loadArchRProject("results")

# amulet
bc_keep = read_tsv(sprintf("%s/post_amulet_barcodes_of_interest_filt.list", FILTER_DIR), col_names=F) %>%
  mutate(X1=sprintf("%s#%s", LIBRARY, X1))
dim(proj@cellColData)

proj = proj[proj$cellNames %in% bc_keep$X1,]

#save(proj, file=sprintf("%s/proj_filt.RData", WORKING_DIR) )

### NOW RUN DOWNSTREAM STEPS ###

proj <- addIterativeLSI(ArchRProj = proj, 
                        useMatrix = "TileMatrix", 
                        name = "IterativeLSI")
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

proj <- addDoubletScores(
  input = proj,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI")

plotEmbedding(ArchRProj = proj, 
              colorBy = "cellColData", 
              name = "Clusters", 
              embedding = "UMAP")

plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      embedding = "UMAP"), name="Clusters")

plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "nFrags", 
                      embedding = "UMAP"), name="nFrags")
save(proj, file=sprintf("%s/proj2.RData", WORKING_DIR))


# markers
proj_markers <- getMarkerFeatures(proj)
save(proj_markers, file=sprintf("%s/proj_markers.RData", WORKING_DIR))

markerList <- getMarkers(proj_markers, 
                         cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
                         n=100)

top5_by_logfc <- lapply(markerList@listData, function(x){
  head(x[order(-x$Log2FC), ], 5)
})

top5_by_fdr <- lapply(markerList@listData, function(x){
  head(x, 5)
})

proj_markers_top5 <- proj_markers[proj_markers@elementMetadata$name %in% 
                                    unique(do.call(rbind, top5_by_logfc)$name), ]
proj_heatmap <- plotMarkerHeatmap(proj_markers_top5, 
                                  #cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
                                  nLabel=10, 
                                  nPrint=10, 
                                  transpose = TRUE)

plotPDF(ComplexHeatmap::draw(proj_heatmap, 
                             heatmap_legend_side = "bot", 
                             annotation_legend_side = "bot"),
        width=8,
        height=6,
        name = sprintf('%s_heatmap_top5_by_logfc.pdf', LIBRARY))

proj_markers_top5 <- proj_markers[proj_markers@elementMetadata$name %in% 
                                    unique(do.call(rbind, top5_by_fdr)$name), ]
proj_heatmap <- plotMarkerHeatmap(proj_markers_top5, 
                                  #cutOff = "FDR <= 0.01 & Log2FC >= 1.5", 
                                  nLabel=10, 
                                  nPrint=10, 
                                  transpose = TRUE)

plotPDF(ComplexHeatmap::draw(proj_heatmap, 
                             heatmap_legend_side = "bot", 
                             annotation_legend_side = "bot"),
        width=8,
        height=6,
        name = sprintf('%s_heatmap_top5_by_fdr.pdf', LIBRARY))

proj <- addImputeWeights(proj)
saveArchRProject(proj)

