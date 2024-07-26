# Input arguments:
# 1. H5 file name
# 2. min number features a valid cell barcode must have.
# 3. min number of cells a valid feature must be expressed in.
library(Seurat)
args = commandArgs(trailingOnly=TRUE)
CELLR_H5_PATH=args[1]
MINFEATURE=as.numeric(args[2])
MINCELL=as.numeric(args[3])
count=Read10X_h5(CELLR_H5_PATH)


print_message <- function(...) {
  cat("[", format(Sys.time()), "] ", ..., "\n", sep="")
}
print_message(sprintf(
  "
  -----------
  Running seurat post-filter with the following params:
    CELLR_H5_PATH=%s
    MINFEATURE=%s
    MINCELL=%s
  -----------
  ", CELLR_H5_PATH, MINFEATURE, MINCELL)
)

# may or may not be in a list with CITE-seq data
if (typeof(count) == "list"){
  sobj=CreateSeuratObject(counts = count$`Gene Expression`, project = "temp", 
                          min.cells = MINCELL, min.features = MINFEATURE)
} else {
  sobj=CreateSeuratObject(counts = count, project = "temp", min.cells = mincell, min.features = minfeat)
}
print(dim(sobj))
bc=row.names(sobj@meta.data)
print(bc)
write.table(bc, file="barcodes_of_interest.filt.list", quote=F, row.names=F, col.names=F)
