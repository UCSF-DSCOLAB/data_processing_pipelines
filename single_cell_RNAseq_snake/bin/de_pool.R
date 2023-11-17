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