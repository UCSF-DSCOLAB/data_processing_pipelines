library(ArchR)
library(tidyverse)

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR="/krummellab/data1/immunox/SICCA1/data/single_nuclear_ATAC/processed"
WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)
FMX_DIR=sprintf("%s/%s/freemuxlet/", PARENT_DIR, LIBRARY)

dir.create(WORKING_DIR)
setwd(WORKING_DIR)
addArchRThreads(threads = 4)
addArchRGenome("hg38")


ArrowFiles <- createArrowFiles(
  inputFiles = sprintf("%s/%s/cellranger/fragments.tsv.gz", PARENT_DIR, LIBRARY),
  sampleNames = LIBRARY,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "results",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)
saveArchRProject(proj)

## subset to the cells of interest
fmx_data = read_tsv(sprintf("%s/%s.clust1.samples.gz", FMX_DIR, LIBRARY)) %>%
  filter(DROPLET.TYPE=="SNG") %>%
  mutate(BARCODE=sprintf("%s#%s", LIBRARY, BARCODE)) %>%
  mutate(fmx_cluster=sprintf("CLUST%s", str_extract(BEST.GUESS, "^[0-9]"))) %>%
  column_to_rownames("BARCODE")

proj = subsetCells(proj, rownames(fmx_data))

proj = addCellColData(ArchRProj = proj, 
                      data = fmx_data$fmx_cluster,
                      cells = rownames(fmx_data),
                      name='fmx_cluster')
saveArchRProject(proj)
## now run next steps
proj <- addIterativeLSI(ArchRProj = proj, 
                        useMatrix = "TileMatrix", 
                        name = "IterativeLSI",
                        force=T)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force=T
)


proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI",
                    force=T)

plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      embedding = "UMAP"),
        name= "Clusters")


plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "nFrags", 
                      embedding = "UMAP"),
        name="nFrags")



saveArchRProject(proj)
# top markers

proj_markers <- getMarkerFeatures(proj)
save(proj_markers, file="proj_markers.RData")

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
        name = sprintf('heatmap_top5_by_logfc.pdf'))

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
        name = sprintf('heatmap_top5_by_fdr.pdf'))





