library(tidyverse)
library(gplots)

args = commandArgs(trailingOnly=T)
POOL = args[1]
gtcheck_out = args[2]

tab_in = read_tsv(gtcheck_out, comment="#", col_names=F) 
if (ncol(tab_in) == 6 & "CN" %in% tab_in$X1){ ### for 1.3.1
  my_tab =  tab_in %>%
    filter(X1=="CN") %>%
    filter(!str_detect(X5, "CLUST"), str_detect(X6, "CLUST")) %>%
    select(-X1, -X3, -X4) %>%
    mutate(X5=str_replace_all(X5, "_[0-9]+", "")) %>%
    dplyr::rename(individual=X5, cluster=X6) %>%
    mutate(err=as.numeric(X2)) %>%
    select(-X2)
} else if (ncol(tab_in)==6 & "DC" %in% tab_in$X1){ # 1.18
  my_tab = tab_in %>% 
    select(-X1) %>%
    filter(X6!=0, str_detect(X3, "CLUST"), !str_detect(X2, "CLUST")) %>%
    rename(individual=X2, cluster=X3, err=X4) %>%
    select(-X5, -X6)
} else { # 1.10
  my_tab = tab_in %>%
    filter(!str_detect(X4, "CLUST"), str_detect(X5, "CLUST"), X1=="ERR") %>%
    select(-X1, -X3) %>%
    mutate(X4=str_replace_all(X4, "_[0-9]+", "")) %>%
    dplyr::rename(individual=X4, cluster=X5) %>%
    mutate(err=as.numeric(X2)) %>%
    select(-X2)
}


my_mat = my_tab %>% pivot_wider(names_from="cluster", values_from="err", values_fill=NA) %>%
  column_to_rownames("individual") %>%
  as.matrix() %>%
  t() 
if(ncol(my_mat)==1 | nrow(my_mat) == 1 | any(is.na(my_mat))){
  print("Warning - unable to plot because too few values in heatmap.")
  
} else {
  pdf(sprintf("%s_gtcheck_out.pdf", POOL), height=5, width=5)
  heatmap.2(my_mat, trace="none", scale="none", cexRow = 0.5, cexCol = 0.5)
  dev.off()
}


cutoff = mean(my_tab$err) - sd(my_tab$err)



filt_tab = my_tab %>% 
  filter(err < cutoff) 
filt_tab %>%
  write_tsv(sprintf("%s_gtcheck_map_to_ind.tsv", POOL))

if (length(unique(filt_tab$individual)) != nrow(filt_tab) | 
    length(unique(filt_tab$cluster)) != nrow(filt_tab) | 
    nrow(filt_tab) != nrow(my_mat) ){
  print("Error, the number of assigned samples and clusters does not match the total")
}