library(ArchR)
library(tidyverse)

PARENT_DIR="/dscolab/data2/AUTOIPI/data/single_nuclear_ATAC/processed"
SCRATCH_DIR = "/krummellab/data1/DSCoLab/AUTOIPI/scratch_data/atac_demux/"

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
POOL = strsplit(LIBRARY, "-")[[1]][3]
WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)
DMX_DIR=sprintf("%s/%s/demuxlet/", PARENT_DIR, LIBRARY)
AMULET_DIR=sprintf("%s/%s/amulet/", PARENT_DIR, LIBRARY)

setwd(WORKING_DIR)
addArchRThreads(threads = 4)
addArchRGenome("hg38")

proj = loadArchRProject("results")

## TODO FILTER
cuts = read_yaml(sprintf("%s/cutoffs.yml", WORKING_DIR))
if (!as.logical(cuts$reviewed)){
  print("Error, qc cutoffs are not reviewed!")
} else {
  pre_filt <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
  post_filt = pre_filt %>% 
    filter(TSSEnrichment > cuts$tss.min,
           `log10.nFrags.` > log10(as.numeric(cuts$nFrag.min)))
  
  write(sprintf("Started with %s nuclei, filtered to %s nuclei based on selected cutoffs.", 
                nrow(pre_filt), nrow(post_filt)), file=sprintf("%s/filter_log.txt", WORKING_DIR))
  
  # now replot
  p = ggPoint(
    x = pre_filt$`log10.nFrags.`,
    y = pre_filt$TSSEnrichment, 
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(pre_filt$`log10.nFrags.`, probs = 0.99)),
    ylim = c(0, quantile(pre_filt$TSSEnrichment, probs = 0.99))) + 
    geom_hline(yintercept = cuts$tss.min, lty = "dashed") + 
    geom_vline(xintercept = log10(as.numeric(cuts$nfrag.min)), lty = "dashed")
  pdf(sprintf("%s/filt_tss_vs_nFrag.pdf", WORKING_DIR), height=5, width=5)
  p %>% ggMarginal(type="density")
  dev.off()
}
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
proj <- addDoubletScores(
  input = proj,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE
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

proj = addImputeWeights(proj)
saveArchRProject(proj)



