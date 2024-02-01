# A few assertions to check de-pooling happens successfully

# Original files

scg1_path = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG1/automated_processing/TEST-POOL-DM1-SCG1_raw.rds"
scg2_path = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG2/automated_processing/TEST-POOL-DM1-SCG2_raw.rds"
scg1_obj <- readRDS(scg1_path)
scg2_obj <- readRDS(scg2_path)

# Read in depooled RDS object

path_to_depooled_obj_1 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1/293T_RTG.rds"
path_to_depooled_obj_2 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1/jurkat.rds"
depooled_obj_0 <- readRDS(path_to_depooled_obj_1)
depooled_obj_1 <- readRDS(path_to_depooled_obj_2)

# Meta-data checks

## Assert singlets are correctly counted
scg1_singlet_count <- sum(scg1_obj@meta.data$DROPLET.TYPE == "SNG")
scg2_singlet_count <- sum(scg2_obj@meta.data$DROPLET.TYPE == "SNG")

depooled_obj_0_singlet_count <- sum(depooled_obj_0@meta.data$DROPLET.TYPE == "SNG")
depooled_obj_1_singlet_count <- sum(depooled_obj_1@meta.data$DROPLET.TYPE == "SNG")

stopifnot(identical(scg1_singlet_count + scg2_singlet_count, depooled_obj_0_singlet_count + depooled_obj_1_singlet_count))

# Raw count checks

## Assert all cell_ids in the meta_data object correspond to a raw_count matrix

depooled_obj_0_meta_cell_ids <- rownames(depooled_obj_0@meta.data)
depooled_obj_0_cell_ids <- colnames(depooled_obj_0@assays$RNA@counts)
stopifnot(identical(depooled_obj_0_meta_cell_ids, depooled_obj_0_cell_ids))

depooled_obj_1_meta_cell_ids <- rownames(depooled_obj_1@meta.data)
depooled_obj_1_cell_ids <- colnames(depooled_obj_1@assays$RNA@counts)
stopifnot(identical(depooled_obj_1_meta_cell_ids, depooled_obj_1_cell_ids))