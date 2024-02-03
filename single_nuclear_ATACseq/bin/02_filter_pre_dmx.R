library(yaml) 
library(tidyverse)
library(ArchR)
set.seed(1)

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR="/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed"
OUT_DIR=sprintf("%s/%s/cell_filter/", PARENT_DIR, LIBRARY)
WORKING_DIR=sprintf("%s/%s/archR/", PARENT_DIR, LIBRARY)

pre_filt = read_csv(sprintf("%s/cell_qc_meta.csv", OUT_DIR))

## read cutoffs
cuts = read_yaml(sprintf("%s/cutoffs.yml", OUT_DIR))

if (!as.logical(cuts$reviewed)){
  print("Error, qc cutoffs are not reviewed!")
} else {
  post_filt = pre_filt %>% 
    filter(TSSEnrichment > cuts$tss.min,
           `log10.nFrags.` > log10(as.numeric(cuts$nFrag.min)))
  
  write(sprintf("Started with %s nuclei, filtered to %s nuclei based on selected cutoffs.", 
                nrow(pre_filt), nrow(post_filt)), file=sprintf("%s/filter_log.txt", OUT_DIR))
  
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
  pdf(sprintf("%s/filt_tss_vs_nFrag.pdf", OUT_DIR), height=5, width=5)
  p
  dev.off()
  
  if (file.exists(sprintf("%s/%s/cellranger/singlecell.csv", PARENT_DIR, LIBRARY))){
    unfilt_input = read_csv(sprintf("%s/%s/cellranger/singlecell.csv", PARENT_DIR, LIBRARY))
    
  } else { # for multiome
    unfilt_input = read_csv(sprintf("%s/%s/cellranger/per_barcode_metrics.csv", PARENT_DIR, library))
    atac_cols = colnames(unfilt_input)[sapply(colnames(unfilt_input), function(x) str_detect(x,"atac"))]
    # there are multiple barcodes for multiome -- we use the general one, not the "ATAC" one
    unfilt_input2 = unfilt_input %>% dplyr::select(barcode, setdiff(atac_cols, "atac_barcode"), is_cell, excluded_reason)
    colnames(unfilt_input2) = str_replace_all(colnames(unfilt_input2) , "atac_", "")
    colnames(unfilt_input2) = str_replace_all(colnames(unfilt_input2) , "_reads", "")
    unfilt_input2.1 = unfilt_input2 %>% dplyr::rename(total=raw, duplicate=dup,
                                                      is__cell_barcode=is_cell) 
    list_cols = c("barcode","total","duplicate", "chimeric","unmapped", "lowmapq", "mitochondrial", "nonprimary","passed_filters",
                  "is__cell_barcode","excluded_reason","TSS_fragments", "DNase_sensitive_region_fragments",
                  "enhancer_region_fragments","promoter_region_fragments","on_target_fragments","blacklist_region_fragments",
                  "peak_region_fragments","peak_region_cutsites")  
    setdiff(list_cols, colnames(unfilt_input2.1))
    
    unfilt_input = unfilt_input2.1 %>% mutate(nonprimary=NA, passed_filters=NA, 
                                              DNase_sensitive_region_fragments=NA,
                                              enhancer_region_fragments=NA,
                                              promoter_region_fragments=NA,
                                              on_target_fragments=NA,
                                              blacklist_region_fragments=NA) %>%
      dplyr::select(list_cols)
  }
  
  filt_barcodes = unfilt_input %>% dplyr::filter(barcode %in% pre_filt$cell_id)
  filt_barcodes %>% write_csv(sprintf("%s/pre_amulet_barcodes_filt.csv", OUT_DIR))    

}


