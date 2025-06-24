library(ArchR)
library(tidyverse)
library(ggExtra)
library(ggplot2)
set.seed(1)

NTRHEADS=4
DEFAULT_TSS_MIN=2
DEFAULT_NFRAG_MIN=100

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
FRAGMENTS=args[2]
AMULET_BC=args[3]
DEMUX_OUT=args[4]

dir.create("archR/", showWarnings=T)
setwd("archR/")
addArchRThreads(threads = NTHREADS) 
addArchRGenome("hg38")

cutoffs = tibble("params"=("tss.min", "nFrag.min", "reviewed"), 
  "vals"= c(DEFAULT_TSS_MIN, DEFAULT_NFRAG_MIN, "FALSE"))
cutoffs %>% write_csv(sprintf("%s_cutoffs.csv", LIBRARY))


ArrowFiles <- createArrowFiles(
  inputFiles = FRAGMENTS,
  sampleNames = LIBRARY,
  minTSS = DEFAULT_TSS_MIN,
  minFrags = DEFAULT_NFRAG_MIN, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "results/",
  copyArrows = TRUE 
)
saveArchRProject(proj)


df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

df %>%
  as_tibble(rownames="BARCODE") %>%
  separate(BARCODE, into=c("LIBRARY", "cell_id"), sep="#", remove=F) %>%
  write_csv("cell_qc_meta.csv")

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
pdf("unfilt_tss_vs_nFrag.pdf", height=5, width=5)
p %>% ggMarginal(type="density")
dev.off()




## read in amulet & dmx data subset to the cells of interest and replot
bc_keep = read_tsv(AMULET_BC, col_names=F) %>%
  mutate(X1=sprintf("%s#%s", LIBRARY, X1))

proj = proj[proj$cellNames %in% bc_keep$X1,]

dmx_data = read_tsv(DEMUX_OUT) %>%
  filter(DROPLET.TYPE=="SNG") %>% 
  separate(BARCODE, into=c("LIBRARY", "cell_id"), sep="#", remove=F) %>%
  column_to_rownames("cell_id")

proj = subsetCells(proj, rownames(dmx_data))

proj = addCellColData(ArchRProj = proj, 
                      data = dmx_data$BEST.GUESS,
                      cells = rownames(dmx_data),
                      name='BEST.GUESS')
saveArchRProject(proj)

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

df %>%
  as_tibble(rownames="BARCODE") %>%
  separate(BARCODE, into=c("LIBRARY", "cell_id"), sep="#", remove=F) %>%
  write_csv("post_demux_cell_qc_meta.csv")

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
pdf("post_demux_tss_vs_nFrag.pdf", height=5, width=5)
p %>% ggMarginal(type="density")
dev.off()



