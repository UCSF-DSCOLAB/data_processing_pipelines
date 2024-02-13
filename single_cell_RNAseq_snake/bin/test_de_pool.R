# A few assertions to check de-pooling happens successfully

# Original files

scg1_path = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG1/automated_processing/TEST-POOL-DM1-SCG1_raw.rds"
scg2_path = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG2/automated_processing/TEST-POOL-DM1-SCG2_raw.rds"
scg1_obj <- readRDS(scg1_path)
scg2_obj <- readRDS(scg2_path)

# Read in depooled RDS object

path_to_depooled_obj_1 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1_SCG1_raw/293T_RTG.rds"
path_to_depooled_obj_2 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1_SCG1_raw/jurkat.rds"
path_to_depooled_obj_3 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1_SCG2_raw/293T_RTG.rds"
path_to_depooled_obj_4 = "/krummellab/data1/amazzara/tutorial_lib_sep/data/single_cell_GEX/processed/depooled/DM1_SCG2_raw/jurkat.rds"

depooled_obj_1 <- readRDS(path_to_depooled_obj_1)
depooled_obj_2 <- readRDS(path_to_depooled_obj_2)
depooled_obj_3 <- readRDS(path_to_depooled_obj_3)
depooled_obj_4 <- readRDS(path_to_depooled_obj_4)

# Meta-data checks

## Assert singlets are correctly counted
scg1_singlet_count <- sum(scg1_obj@meta.data$DROPLET.TYPE == "SNG")
scg2_singlet_count <- sum(scg2_obj@meta.data$DROPLET.TYPE == "SNG")

depooled_obj_1_singlet_count <- sum(depooled_obj_1@meta.data$DROPLET.TYPE == "SNG")
depooled_obj_2_singlet_count <- sum(depooled_obj_2@meta.data$DROPLET.TYPE == "SNG")
depooled_obj_3_singlet_count <- sum(depooled_obj_3@meta.data$DROPLET.TYPE == "SNG")
depooled_obj_4_singlet_count <- sum(depooled_obj_4@meta.data$DROPLET.TYPE == "SNG")

stopifnot(identical(scg1_singlet_count, depooled_obj_1_singlet_count + depooled_obj_2_singlet_count))
stopifnot(identical(scg2_singlet_count, depooled_obj_3_singlet_count + depooled_obj_4_singlet_count))

# Raw count checks

## Assert all cell_ids in the meta_data object correspond to a raw_count matrix

depooled_obj_1_meta_cell_ids <- rownames(depooled_obj_1@meta.data)
depooled_obj_1_cell_ids <- colnames(depooled_obj_1@assays$RNA@counts)
stopifnot(identical(depooled_obj_1_meta_cell_ids, depooled_obj_1_cell_ids))

depooled_obj_2_meta_cell_ids <- rownames(depooled_obj_2@meta.data)
depooled_obj_2_cell_ids <- colnames(depooled_obj_2@assays$RNA@counts)
stopifnot(identical(depooled_obj_2_meta_cell_ids, depooled_obj_2_cell_ids))

depooled_obj_3_meta_cell_ids <- rownames(depooled_obj_3@meta.data)
depooled_obj_3_cell_ids <- colnames(depooled_obj_3@assays$RNA@counts)
stopifnot(identical(depooled_obj_3_meta_cell_ids, depooled_obj_3_cell_ids))

depooled_obj_4_meta_cell_ids <- rownames(depooled_obj_4@meta.data)
depooled_obj_4_cell_ids <- colnames(depooled_obj_4@assays$RNA@counts)
stopifnot(identical(depooled_obj_4_meta_cell_ids, depooled_obj_4_cell_ids))