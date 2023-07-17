library(Seurat)
library(tidyverse)
library(ggpubr)
library(dsb)
library(ggExtra)

#? TODO droplet filtering / dualscatter?

args = commandArgs(trailingOnly=T)

LIBRARY=args[1] 
BASE_DIR=args[2]
source(sprintf("%s/bin/seurat_utils.R", BASE_DIR))

sobj = readRDS(sprintf("%s_raw.rds", LIBRARY))
param_file = sprintf("%s_cutoffs.csv", LIBRARY)
params_df = read_csv(param_file)
reviewed_bool = params_df %>% filter(parameter=="reviewed") %>% mutate(value=as.logical(toupper(value))) %>% pull(value)
if (!reviewed_bool){
  message("The cutoffs in the parameter file (", param_file, ") are not reviewed. Please review the cutoffs, set the value of 'reviewed' to 'TRUE' and rerun the pipeline.")
  quit(status=10)
}

adt.present = ("ADT" %in% names(sobj@assays))

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

# plots with final reviewed params
params = load_params(params_df) 
plot_list = make_plots(sobj, params)
num_rows = ifelse(adt.present, 3, 2)
pdf( sprintf("%s_diagnostic_plots_final.pdf", LIBRARY) , width=25, height=7*num_rows)
ggarrange(plotlist=plot_list, ncol=3,nrow=num_rows)
dev.off()


# Adding Seurat parameters to the misc slot of the Seurat object.
sobj@misc$scStat$seurat_params = params

# Filter out the low-quality cells.
sobj = filter_cells(sobj, params, adt.present)

# Adding cell count to misc slot after filtering low-quality cells.
sobj@misc$scStat$nCells_after_filter = ncol(sobj)

list_cells = colnames(sobj)
tibble("cell"=list_cells) %>% write_tsv( file="barcodes_of_interest.filt.list", col_names=FALSE)

### now run the next steps of process_10x_with_seurat

