library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)
library(gridExtra)
library(cowplot)
library(dittoSeq)

args = commandArgs(trailingOnly=T)

LIBRARY=args[1] 
BASE_DIR=args[2]
CELLR_H5_PATH=args[3]
KEEP_FMX_SNG=(args[4]=="true")
KEEP_ONLY_SNG=(args[5]=="true")
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

print_message(sprintf(
  "
  -----------
  Running seurat post-filter with the following params:
    LIBRARY=%s
    BASE_DIR=%s
    CELLR_H5_PATH=%s
    KEEP_FMX_SNG=%s
    KEEP_ONLY_SNG=%s
  -----------
  ", LIBRARY, BASE_DIR, CELLR_H5_PATH, KEEP_FMX_SNG, KEEP_ONLY_SNG)
)


# params: TODO make these parameters
VARS_TO_REGRESS=c("percent.mt","percent.ribo","nCount_RNA","nFeature_RNA",
  "S.Score", "G2M.Score")
truncate.adt = TRUE
sct.transform = FALSE


sobj = readRDS(sprintf("%s_raw.rds", LIBRARY))
param_file = sprintf("%s_cutoffs.csv", LIBRARY)
params_df = read_csv(param_file)
reviewed_bool = params_df %>% filter(parameter=="reviewed") %>% mutate(value=as.logical(toupper(value))) %>% pull(value)
if (!reviewed_bool){
  message("The cutoffs in the parameter file (", param_file, ") are not reviewed. Please review the cutoffs, set the value of 'reviewed' to 'TRUE' and rerun the pipeline.")
  quit(status=10)
}

adt.present = ("ADT" %in% names(sobj@assays))

# plots with final reviewed params
params = load_params(params_df)

plot_list = suppressWarnings(make_plots(sobj, params, adt.present))
num_rows = ifelse(adt.present, 5, 3)
merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
ggsave(merge, file=sprintf("%s_diagnostic_plots_reviewed.png", LIBRARY) , width=30, height=7*num_rows, bg="white", dpi=72)


# Adding Seurat parameters to the misc slot of the Seurat object.

sobj@misc$scStat$seurat_params = params
sobj@misc$scStat$nCells_before_filter = ncol(sobj)

print(sprintf("Before filtering there are %s cells", ncol(sobj)))


# Filter out the low-quality cells.
if ("DROPLET.TYPE.FINAL" %in% colnames(sobj@meta.data)){
  no_nas = subset(sobj, DROPLET.TYPE.FINAL %in% c("AMB", "Intra.DBL", "Inter.DBL", "SNG"))
  plot_list = suppressWarnings(make_plots(no_nas, params, adt.present, group="DROPLET.TYPE.FINAL"))
} else {
  no_nas = subset(sobj, DROPLET.TYPE %in% c("AMB", "DBL", "SNG"))
  plot_list = suppressWarnings(make_plots(no_nas, params, adt.present, group="DROPLET.TYPE"))
}
num_rows = ifelse(adt.present, 5, 3)
merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
ggsave(sprintf("%s_doublet_plot_reviewed.png",LIBRARY),  width=30, height=7*num_rows, bg="white", dpi=72)

sobj = filter_cells(sobj, params, adt.present)
if (KEEP_FMX_SNG){
  sobj = subset(sobj, DROPLET.TYPE=="SNG")
}
if (KEEP_ONLY_SNG){
  if ("DROPLET.TYPE.FINAL" %in% colnames(sobj@meta.data)){
    sobj = subset(sobj, DROPLET.TYPE.FINAL=="SNG")
  } else {
    sobj = subset(sobj, DROPLET.TYPE=="SNG")
  }
}
print(sprintf("Following filtering there are %s cells remaining", ncol(sobj)))


df = quantile_frac_table(sobj, adt.present)
write.table(df, sprintf("%s_quantiles_post.tsv", LIBRARY), sep="\t")

# final plots 
plot_list = suppressWarnings(make_plots(sobj, params, adt.present, add_stats=F))
num_rows = ifelse(adt.present, 5, 3)
merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
ggsave(merge, file=sprintf("%s_diagnostic_plots_filtered.png", LIBRARY) , width=30, height=7*num_rows, bg="white", dpi=72)


# Adding cell count to misc slot after filtering low-quality cells.
sobj@misc$scStat$nCells_after_filter = ncol(sobj)

# Preparing cell-cycle genes. cc.genes is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


# Normalize data for cell-cycle scoring
sobj = sobj %>% NormalizeData()


sobj = CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DefaultAssay(sobj) = "RNA"

saveRDS(sobj, file=sprintf("%s_filtered.rds", LIBRARY))

# write out the list of filtered barcodes


# Normalize and Scale the data while regressing out Mito content, Ribo content, number of genes, UMIs, and cell cycle scores.

# SCTransform if desired
if (sct.transform){
  sobj = sobj %>% 
    SCTransform(vars.to.regress = VARS_TO_REGRESS,
    return.only.var.genes=FALSE,
    verbose=FALSE)
} else {
  sobj = sobj %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = VARS_TO_REGRESS )
}


# normalize ADT data with DSB if its present
if (adt.present){
   raw = Read10X_h5(CELLR_H5_PATH)
   #background_drops <- colnames(raw$`Gene Expression`)[colSums(raw$`Gene Expression`)<100]
   
   # update background calculation
   md = data.frame(
    prot.size = log10(colSums(raw$`Antibody Capture`)), # -Inf from 0 ones
    rna.size = log10(colSums(raw$`Gene Expression`)) # -Inf from 0 ones
   )   
   md = md[md$rna.size > 0 & md$prot.size > 0, ]
   if (!any(str_detect(rownames(params), "background"))){
    background_drops = rownames(
      md[ md$prot.size > 1.5 & 
      md$rna.size < max(2.5, log10(params['nFeature_RNA.lower',]+0.1)), ]
    )
   } else {
    background_drops = rownames(
      md[ md$prot.size > log10(params['background_ADT.lower',]) & 
        md$prot.size < log10(params['background_ADT.upper',]) &
      md$rna.size < log10(params['background_RNA.upper',]), ]
    )
   }
    
   print(sprintf("FILTER background with %s empty drops",
                  length(background_drops)))
   
   ADT_background = raw$`Antibody Capture`[, background_drops]
   adt_norm = process_adt(sobj[['ADT']]@counts, ADT_background)

   # truncate values below zero 
   if (truncate.adt){
     adt_norm[adt_norm < 0] = 0
   }
   sobj@assays$ADT@data = adt_norm

   # additional ADT diagnostic plots
   adt_plot(sobj)
   ggsave(sprintf("%s_normalized_ADT_dist.png", LIBRARY), height=30, width=30)
}

sobj = sobj %>%
  RunPCA(verbose=F) %>%
  RunUMAP(dims=1:30, verbose=F)
sobj = sobj %>%
  FindNeighbors(dims=1:30, k.param=20, verbose=F) 

sobj = sobj %>%
  FindClusters(verbose=F, algorithm=4, method="igraph")

DimPlot(sobj, group.by="seurat_clusters")
ggsave(sprintf("%s_umap.pdf", LIBRARY), height=5, width=5)

# confound PCA, UMAP
FeaturePlot(sobj, features=VARS_TO_REGRESS, reduction="pca")
ggsave(sprintf("%s_confound_pca.pdf", LIBRARY), height=9, width=6)

FeaturePlot(sobj, features=VARS_TO_REGRESS, reduction="umap")
ggsave(sprintf("%s_confound_umap.pdf", LIBRARY), height=9, width=6)


# SAMPLE.By.SNPs
sobj@meta.data$SAMPLE.by.SNPs = as.factor(sapply(as.vector(sobj@meta.data$BEST.GUESS), function(x) {
  temp <- paste(sort(unique(strsplit(x, ',')[[1]])), collapse='_')
}))

DimPlot(sobj, group.by="SAMPLE.by.SNPs")
ggsave(sprintf("%s_samples_by_SNPs_umap.pdf", LIBRARY), height=5, width=5)

FeaturePlot(sobj, "doubletFinderScore")
ggsave(sprintf("%s_doubletScore_umap.pdf", LIBRARY), height=5, width=5)

# immune triplot
if ("PTPRC" %in% rownames(sobj)){
  pdf(sprintf("%s_immune_triplot_umap.pdf", LIBRARY), height=6, width=7)
  TriPlot(sobj, features="PTPRC")
  dev.off()
}
# Save the final seurat object
saveRDS(sobj, file=sprintf("%s_processed.rds", LIBRARY))

list_cells = colnames(sobj)
tibble("cell"=list_cells) %>% write_tsv( file="barcodes_of_interest_filt.list", col_names=FALSE)


