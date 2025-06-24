library(tidyverse)
args = commandArgs(trailingOnly=T)
bc_list=args[1]

mult_bc = read_tsv("amulet/MultipletBarcodes_01.txt"), col_names=F) %>%
  dplyr::rename(cell_id=X1) %>%
  mutate(amuletDBL=TRUE)

pre_filt = read_csv(bc_list)
post_filt = pre_filt %>%
  filter(!barcode %in% mult_bc$cell_id)
write(sprintf("Removing %s amuletDBLs, resulting in %s nuclei.",
              nrow(mult_bc), nrow(post_filt)), file="amulet/filter_log.txt",
      append=T)

write_tsv(post_filt %>% 
            dplyr::select(barcode), col_names=F, 
          file="amulet/post_amulet_barcodes_of_interest.filt.list")