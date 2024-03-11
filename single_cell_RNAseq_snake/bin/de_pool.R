library(Seurat)
depool_and_save <- function(input_file, output_dir) {
  if (!file.exists(input_file)) {
    stop("File not found: ", input_file)
  }

  # Load the Seurat object
  seurat_obj <- readRDS(input_file)

  # Determine sample ids
  seurat_obj$sample_id <- sapply(strsplit(as.character(seurat_obj$BEST.GUESS), ","), function(x) x[1])

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Identify unique sample IDs and save subsets
  unique_ids <- unique(seurat_obj$sample_id)
  for (id in unique_ids) {
    cells_sng <- WhichCells(seurat_obj, expression = sample_id == id & seurat_obj$DROPLET.TYPE == "SNG")
    if (length(cells_sng) > 0) {
      subset_obj <- subset(seurat_obj, cells = cells_sng)
      filename <- file.path(output_dir, sprintf("%s.rds", id))
      saveRDS(subset_obj, file = filename)
    } else{
        warning("No singlet cells exist for: ", id, " in input_file ", input_file)
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

cat("Input file:", input_file, "\n")
cat("Output dir:", output_dir, "\n")

cat("Depooling seurat_objects...\n" )
depool_and_save(input_file, output_dir)
