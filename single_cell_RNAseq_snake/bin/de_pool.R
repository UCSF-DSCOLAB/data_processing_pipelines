# This script depools demultiplexed single cell RNA seq data

load_seurat_objects <- function(file_paths) {
    seurat_objects <- list()

    for (file_path in file_paths) {
        if (file.exists(file_path)) {
            seurat_objects[[basename(file_path)]] <- readRDS(file_path)
        } else {
            stop("File not found: ", file_path)
        }
    }
    return(seurat_objects)
}


save_seurat_objects <- function(output_dir, seurat_list) {
  # Check if the provided list is actually a list of Seurat objects
  if (!is.list(seurat_list)) {
    stop("The provided 'seurat_list' is not a list.")
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
  }

  # Iterate through each Seurat object in the list
  for (id in names(seurat_list)) {
    obj <- seurat_list[[id]]

    # Define the filename using the output directory and the ID
    filename <- file.path(output_dir, paste0(id, ".rds"))

    # Save the Seurat object
    saveRDS(obj, file = filename)
  }
}

depool_by_id <- function(seurat_objects) {
  depooled_objs <- list()
  
  for (seurat_obj in seurat_objects) {
    # Extract metadata
    metadata <- seurat_obj@meta.data
    
    # Split BEST.GUESS values and extract first for SNG, then store in 'sample_id'
    metadata$sample_id <- sapply(strsplit(as.character(metadata$BEST.GUESS), ","), function(x) x[1])

    # Update the metadata in the Seurat object
    seurat_obj <- Seurat::AddMetaData(object = seurat_obj, metadata = metadata)

    # Create or update subsets for each unique value in 'sample_id'
    unique_ids <- unique(metadata$sample_id)
    
    for (id in unique_ids) {
      # Select only singlets cells for this sample_id
      cells_sng <- rownames(metadata[metadata$sample_id == id & metadata$DROPLET.TYPE == "SNG", ])
      
      # Proceed only if there are singlet cells for this sample_id
      if (length(cells_sng) > 0) {
        if (id %in% names(depooled_objs)) {
          # Merge if subset already exists
          depooled_objs[[id]] <- merge(depooled_objs[[id]], subset(seurat_obj, cells = cells_sng))
        } else {
          # Create new subset if it doesn't exist
          depooled_objs[[id]] <- subset(seurat_obj, cells = cells_sng)
        }
      }
    }
  }
  return(depooled_objs)
}


args <- commandArgs(trailingOnly = TRUE)

# All arguments except the last one are input files
input_files <- args[1:(length(args) - 1)]
# The last argument is the output file
output_dir <- args[length(args)]

cat("Input files:", input_files, "\n")
cat("Output dir:", output_dir, "\n")
cat("Loading seurat objects... \n")
targets <- load_seurat_objects(input_files)
cat("Depooling seurat_objects...\n" )
depooled_objects <- depool_by_id(targets)
cat("Saving depooled objects...\n")
save_seurat_objects(output_dir, depooled_objects)