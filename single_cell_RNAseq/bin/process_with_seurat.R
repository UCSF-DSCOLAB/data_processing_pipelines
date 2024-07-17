library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)

args = commandArgs(trailingOnly=TRUE)
LIBRARY = args[1]
main_data_type = args[2] # CITE or GEX
SOBJ = args[3]
BASE_DIR=args[4]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))


# Input arguments:
# 1. LIBRARY (e.g. SICCA1-POOL-GN1-SCG1)
#     getting the sobj from doubletfinder
# 2. data_types (SCG, CITE, TCR, BCR)
# 3. where to look for the TCR/BCR data -- and has it completed?

### PARAMS
# default QC cutoffs
# mitochondia gene name pattern
# ribosomal gene name pattern
QCCUTOFFS = args[5]
MTPATTERN = "^MT-"
RIBOPATTERN = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"

CELLR_H5_PATH=args[6]



# Process the qc cutoff values.
params_df = read_csv(QCCUTOFFS)
if (!file.exists( sprintf("%s_cutoffs.csv", LIBRARY))){
  write_csv(params_df, sprintf("%s_cutoffs.csv", LIBRARY))
}

# Read in the Seurat object
sobj = readRDS(SOBJ)

# check if we have ADT data too
adt.present = (main_data_type=="CITE")

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
    rm(ab_data, adt_assay)
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

  if (length(isotype_ctls) > 0){
    print("the following will be used as isotype controls:")
    print(paste(isotype_ctls, collapse=","))
    isotype_ctl_data = as.matrix(ADT_counts[isotype_ctls,])
    sobj = AddMetaData(sobj, apply(isotype_ctl_data, 2, max), "isotype_ctl_max")
  } else {
    print("Warning - no isotype controls found in ADT data")
    sobj = AddMetaData(sobj, rep(0, ncol(sobj)), "isotype_ctl_max")
  }

}

saveRDS_(sobj, file=sprintf("%s_raw.rds", LIBRARY))

# Generate various diagnostic plots
params = load_params(params_df) 

# make tsvs detailing the quantile fractions
df = quantile_frac_table(sobj, adt.present)
write.table(df, sprintf("%s_quantiles_pre.tsv", LIBRARY), sep="\t")


plot_list = make_plots(sobj, params)
num_rows = ifelse(adt.present, 3, 2)
pdf( sprintf("%s_diagnostic_plots_pre.pdf", LIBRARY) , width=25, height=7*num_rows)
ggarrange(plotlist=plot_list, ncol=3,nrow=num_rows)
dev.off()
