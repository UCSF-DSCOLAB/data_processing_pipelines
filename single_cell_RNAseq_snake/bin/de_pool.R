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


save_seurat_objects <- function(filename, seurat_list) {
  # Check if the provided list is actually a list of Seurat objects
  if (!is.list(seurat_list)) {
    stop("The provided 'seurat_list' is not a list.")
  }

  for (obj in seurat_list) {
    if (!inherits(obj, "Seurat")) {
      stop("All elements of 'seurat_list' must be Seurat objects.")
    }
  }

  dir_to_create <- dirname(filename) # Extract the directory part of the filename

  # Create the directory if it doesn't exist
  if (!dir.exists(dir_to_create)) {
      dir.create(dir_to_create, recursive = TRUE)
  }

  # Save the list of Seurat objects to a file
  saveRDS(seurat_list, file = filename)
}


depool_by_id <- function(seurat_objects) {
  depooled_objs <- list()
  
  for (seurat_obj in seurat_objects) {
    # Extract metadata
    metadata <- seurat_obj@meta.data
    
    # Split BEST.GUESS values and extract first for SNG, then store in 'sample_id'
    metadata$sample_id <- sapply(strsplit(as.character(metadata$BEST.GUESS), ","), function(x) x[1])
    seurat_obj@meta.data <- metadata  # Update the metadata in the Seurat object
    
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
output_file <- args[length(args)]

cat("Input files:", input_files, "\n")
cat("Output file:", output_file, "\n")
cat("Loading seurat objects... \n")
targets <- load_seurat_objects(input_files)
cat("Depooling seurat_objects...\n" )
depooled_objects <- depool_by_id(targets)
cat("Saving depooled objects...\n")
save_seurat_objects(output_file, depooled_objects)