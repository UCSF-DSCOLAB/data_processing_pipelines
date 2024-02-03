library(yaml) 
library(tidyverse)
library(ArchR)
set.seed(1)

args = commandArgs(trailingOnly=T)
LIBRARY=args[1]
PARENT_DIR="/krummellab/data1/immunox/AUTOIPI/data/single_nuclear_ATAC/processed"
OUT_DIR=sprintf("%s/%s/cell_filter/", PARENT_DIR, LIBRARY)
AMULET_DIR=sprintf("%s/%s/amulet/", PARENT_DIR, LIBRARY)


mult_bc = read_tsv(sprintf("%s/MultipletBarcodes_01.txt", AMULET_DIR), col_names=F) %>%
  dplyr::rename(cell_id=X1) %>%
  mutate(amuletDBL=TRUE)

pre_filt = read_csv(sprintf("%s/pre_amulet_barcodes_filt.csv", OUT_DIR))
post_filt = pre_filt %>%
  filter(!barcode %in% mult_bc$cell_id)
write(sprintf("Removing %s amuletDBLs, resulting in %s nuclei.",
              nrow(mult_bc), nrow(post_filt)), file=sprintf("%s/filter_log.txt", OUT_DIR),
      append=T)

write_tsv(post_filt %>% 
            dplyr::select(barcode), col_names=F, 
          file=sprintf("%s/post_amulet_barcodes_of_interest_filt.list", OUT_DIR))
