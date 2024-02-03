library(ArchR)
library(tidyverse)
set.seed(1)

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR="/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed"
WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)
AMULET_DIR=sprintf("%s/%s/amulet/", PARENT_DIR, LIBRARY)
DMX_DIR=sprintf("%s/%s/demuxlet/", PARENT_DIR, LIBRARY)

# load
load(file=sprintf("%s/proj.RData", WORKING_DIR))
loadArchRProject("results")

# amulet
mult_bc = read_tsv(sprintf("%s/MultipletBarcodes_01.txt", AMULET_DIR), col_names=F) %>%
  dplyr::rename(cell_id=X1) %>%
  mutate(amuletDBL=TRUE)


# demuxlet
dmx_out = read_tsv(sprintf("%s/%s.clust1.samples.reduced.tsv", DMX_DIR, LIBRARY)) %>%
  column_to_rownames("BARCODE")

proj = proj[dmx_out$BARCODE,] # TODO: check on this filtering
 
proj = addCellColData(ArchRProj = proj, 
                       data = dmx_out$DROPLET.TYPE,
                       cells = rownames(dmx_out),
                       name='DROPLET.TYPE')
proj = addCellColData(ArchRProj = proj, 
                       data = dmx_out$BEST.GUESS,
                       cells = rownames(dmx_out),
                       name='BEST.GUESS')

dmx_counts = table(proj$DROPLET.TYPE)
proj = proj[proj$DROPLET.TYPE=="SNG", ]
write(sprintf("Removing %s dmx DBL and %s dmx AMB, resulting in %s nuclei.",
              dmx_counts[["DBL"]], dmx_counts[["AMB"]], nrow(proj@cellColData)), file=sprintf("filter_log.txt", OUT_DIR),
      append=T)
saveArchRProject(proj)
save(proj, file=sprintf("proj_filt.RData", WORKING_DIR) )

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

plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      embedding = "UMAP"), name="Clusters")

plotPDF(plotEmbedding(ArchRProj = proj, 
                      colorBy = "cellColData", 
                      name = "nFrags", 
                      embedding = "UMAP"), name="nFrags")
saveArchRProject(proj)
save(proj, file=sprintf("%s/proj2.RData", WORKING_DIR))


# markers
proj_markers <- getMarkerFeatures(proj)
save(proj_markers, file=sprintf("%s/proj2_markers.RData", WORKING_DIR))

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
save(proj, file=sprintf("%s/proj2.RData", WORKING_DIR))

