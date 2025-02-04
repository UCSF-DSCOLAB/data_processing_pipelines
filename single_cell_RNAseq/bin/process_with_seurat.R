library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)
library(dittoSeq)

args = commandArgs(trailingOnly=TRUE)
LIBRARY = args[1]
MAIN_DATA_TYPE = args[2] # CITE or GEX
SOBJ = args[3]
BASE_DIR=args[4]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

QCCUTOFFS = args[5] # default QC cutoffs
CELLR_H5_PATH=args[6] # used for processing ADT
CR_FILT_BC=args[7] # used for plotting cellranger barcodes

# Input arguments:
# 1. LIBRARY (e.g. SICCA1-POOL-GN1-SCG1)
#     getting the sobj from doubletfinder
# 2. data_types (SCG, CITE, TCR, BCR)
# 3. where to look for the TCR/BCR data -- and has it completed?

### PARAMS
# mitochondia gene name pattern
# ribosomal gene name pattern
MTPATTERN = "^MT-"
RIBOPATTERN = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"


print_message(sprintf(
  "
  -----------
  Running seurat processing with the following params:
    LIBRARY=%s
    MAIN_DATA_TYPE=%s
    SOBJ=%s
    BASE_DIR=%s
    QCCUTOFFS=%s
    CELLR_H5_PATH=%s
    CR_FILT_BC=%s
  -----------
  ", LIBRARY, MAIN_DATA_TYPE, SOBJ, BASE_DIR, QCCUTOFFS, CELLR_H5_PATH, CR_FILT_BC)
)



# Process the qc cutoff values.
params_df = read_csv(QCCUTOFFS)
params = load_params(params_df) 

if (!file.exists( sprintf("%s_cutoffs.csv", LIBRARY))){
  write_csv(params_df, sprintf("%s_cutoffs.csv", LIBRARY))
}

# Read in the Seurat object
sobj = readRDS(SOBJ)

# check if we have ADT data too
adt.present = (MAIN_DATA_TYPE=="CITE")


if (adt.present){
  # then load it
  raw = Read10X_h5(CELLR_H5_PATH)
  if (typeof(raw) != "list" | ! "Antibody Capture" %in% names(raw)){
    print("No ADT data present. skipping this step")
    rm(raw)
    gc()
    adt.present=FALSE
  } else{
    ab_data=raw$`Antibody Capture`
    # add the ADT counts so we can add the nCount_ADT filter
    adt_assay = CreateAssayObject(counts=ab_data)
    adt_assay = subset(adt_assay, cells=colnames(sobj))
    stopifnot(all(colnames(sobj)==colnames(adt_assay)))
    sobj[["ADT"]] = adt_assay
    #rm(ab_data, adt_assay)
    adt.present=TRUE
  }
} else {
  adt.present=FALSE
}
if (adt.present){
  print(sprintf("ADT data loaded"))
}


# Adding cell count to misc slot before filtering low-quality cells.
sobj@misc$scStat$nCells_before_filter = ncol(sobj)

# Add percent mito and percent ribo to meta.data.
sobj = sobj %>%
  PercentageFeatureSet(pattern = MTPATTERN, col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = RIBOPATTERN, col.name = "percent.ribo")

# add a field for isotype control max
if (adt.present){
  ADT_counts = sobj@assays[["ADT"]]@counts
  isotype_ctls = rownames(ADT_counts)[str_detect(tolower(rownames(ADT_counts)), "iso")]

  plt = background_plot(raw$`Gene Expression`, ab_data, params)

  if (length(isotype_ctls) > 0){
    print("the following will be used as isotype controls:")
    print(paste(isotype_ctls, collapse=","))
    isotype_ctl_data = as.matrix(ADT_counts[isotype_ctls,])
    sobj = AddMetaData(sobj, apply(isotype_ctl_data, 2, max), "isotype_ctl_max")

    # TODO: isotype control plot
    plt_iso = isotype_ctl_plot(ADT_counts, isotype_ctl_data, params)

    ggarrange(plt, plt_iso, ncol=2)
    ggsave(sprintf("%s_adt_diagnostic_plots.png", LIBRARY), height=7, width=14, dpi=72, bg="white")

  } else {
    print("Warning - no isotype controls found in ADT data")
    sobj = AddMetaData(sobj, rep(0, ncol(sobj)), "isotype_ctl_max")

    ggarrange(plt)
    ggsave(sprintf("%s_adt_diagnostic_plots.png", LIBRARY), height=7, width=7, dpi=72, bg="white")

  }

}

saveRDS(sobj, file=sprintf("%s_raw.rds", LIBRARY))
# Generate various diagnostic plots


# make tsvs detailing the quantile fractions
df = quantile_frac_table(sobj, adt.present)
write.table(df, sprintf("%s_quantiles_pre.tsv", LIBRARY), sep="\t")

print("quantiles created and saved")

plot_list = suppressWarnings(make_plots(sobj, params, adt.present))
num_rows = ifelse(adt.present, 5, 3)
merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
ggsave(merge, file=sprintf("%s_diagnostic_plots_default.png", LIBRARY) , width=30, height=7*num_rows, bg="white", dpi=72)



cellranger_bc = read_tsv(CR_FILT_BC, col_names=F) %>% pull(X1)
sobj@meta.data$cellranger_cell = sobj@meta.data %>% 
    as_tibble(rownames="cell_id") %>%
    mutate(cellranger_cell=(cell_id %in% cellranger_bc)) %>%
    pull(cellranger_cell)
  

plot_list = suppressWarnings(make_plots(sobj, params, adt.present, group="cellranger_cell"))
num_rows = ifelse(adt.present, 5, 3)
merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
ggsave(merge, file=sprintf("%s_cr_diagnostic_plots_default.png", LIBRARY) , width=30, height=7*num_rows, bg="white", dpi=72)

# if doublet info present, then plot it
if ("DROPLET.TYPE" %in% colnames(sobj@meta.data)){
  if ("DROPLET.TYPE.FINAL" %in% colnames(sobj@meta.data)){
    no_nas = subset(sobj, DROPLET.TYPE.FINAL %in% c("AMB", "Intra.DBL", "Inter.DBL", "SNG"))
    plot_list = suppressWarnings(make_plots(no_nas, params, adt.present, group="DROPLET.TYPE.FINAL"))
  } else {
    no_nas = subset(sobj, DROPLET.TYPE %in% c("AMB", "DBL", "SNG"))
    plot_list = suppressWarnings(make_plots(no_nas, params, adt.present, group="DROPLET.TYPE"))
  }
  num_rows = ifelse(adt.present, 5, 3)
  merge = ggarrange(plotlist=plot_list, ncol=4,nrow=num_rows)
  ggsave(sprintf("%s_doublet_plot_default.png",LIBRARY),  width=30, height=7*num_rows, bg="white", dpi=72)
}
