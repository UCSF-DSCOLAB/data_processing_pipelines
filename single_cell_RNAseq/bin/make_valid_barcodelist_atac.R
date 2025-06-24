library(tidyverse)

args = commandArgs(trailingOnly=T)

BC_FILE=args[1]
MIN_FRAGS=as.numeric(args[2])
unfilt_input = read_csv(BC_FILE)
filt_barcodes = unfilt_input %>% filter(passed_filters > MIN_FRAGS) %>% 
      filter(barcode != "NO_BARCODE") 
filt_barcodes %>% select(barcode) %>%
    write_tsv("demuxlet/barcodes_of_interest.filt.list", col_names=F)
filt_barcodes %>% write_csv("amulet/pre_amulet_barcodes_filt.csv")


