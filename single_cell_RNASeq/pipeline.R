suppressPackageStartupMessages({
  library(assertthat)
  library(Seurat)
  library(dplyr)    # inline modification of matrices
  library(cowplot)  # Pretty plots
  library(ggplot2)  # Pretty plots
  library(scales)   # to access break formatting functions
  library(yaml)     # For reading cutoffs
})

# Set a seed so results are reproducible
seed = 21212

defaultArgs = list(
  SAMPLE_YML=list(
    value=NA,
    class=as.character,
    description="A file containing the samples to process.",
    global=TRUE
  ),
  WORKING_FOLDER=list(
    value=Sys.getenv("TMPDIR"),
    class=as.character,
    description="A working directory. Run warnings are written here",
    global=TRUE
  ),
  RSCRIPTS_DIR=list(
    value=paste(Sys.getenv("SCRIPTS"), "R", sep="/"),
    class=as.character,
    description="Where to find aux scripts",
    global=TRUE
  ),
  GENESET_DIR=list(
    value=paste(Sys.getenv("KLAB"), "ipi/data/refs/10x/genesets", sep="/"),
    class=as.character,
    description="Where to find genesets for ribo and cell cycling.",
    global=TRUE
  ),
  OUT_FOLDER=list(
    value=NA,
    class=as.character,
    description="Where to write results.",
    global=FALSE
  ),
  SPECIES=list(
    value="human",
    class=as.character,
    description="Species to use for cell cycling and ribo. Must be `human` or `mouse` (case-sensitive).",
    global=FALSE
  ),
  RERUN_STAGE=list(
    value=FALSE,
    class=as.logical,
    description="Should we rerun the last run stage?",
    global=FALSE
  ),
  SECONDARY_ASSAY_NAME=list(
    value="IDX",
    class=as.character,
    description="If there is a secondary assay in the GEX, what should it be called?",
    global=FALSE
  ),
  IDX_ASSAY_NAME=list(
    value="IDX",
    class=as.character,
    description="What should the demultiplexing-related antibody capture assay (if any) be called?",
    global=FALSE
  ),
  IDX_NORMALIZATION_METHOD=list(
    value="CLR",
    class=as.character,
    description="Normalization method to use for the `IDX_ASSAY_NAME` assay",
    global=FALSE
  ),
  ADT_ASSAY_NAME=list(
    value="ADT",
    class=as.character,
    description="What should the other antibody capture assay (if any) be called? Currently this is meant for CITE-Seq",
    global=FALSE
  ),
  ADT_NORMALIZATION_METHOD=list(
    value="CLR",
    class=as.character,
    description="Normalization method to use for the `ADT_ASSAY_NAME` assay",
    global=FALSE
  )
)

# convenience
global_args <- names(defaultArgs)[sapply(defaultArgs, function(x){x$global})]
local_args <- names(defaultArgs)[!names(defaultArgs) %in% global_args]

allowedSampleLevelArgs <- list(
  'GEX_datadir'=as.character,
  'cells_loaded'=as.integer,
  'ADT_datadir'=as.character,
  'IDX_datadir'=as.character,
  'BCR_datadir'=as.character,
  'TCR_datadir'=as.character,
  'multiplet_rate'=as.numeric,
  'metadata'=as.list
  )

print_help <- function(defaultArgs) {
  helptext = "
This is an R Script to process a 10X samples from raw feature/barcode matrices.

USAGE:
  
  Rscript process_10x_with_Seurat.R SAMPLE_YML=/path/to/samples.yml \
                                    OPTIONAL_ARG1=OPTIONAL_VAL1 \
                                    OPTIONAL_ARG2=OPTIONAL_VAL2 \
                                              ....
                                    OPTIONAL_ARGn=OPTIONAL_VALn \

Required Arguments:
SAMPLE_YML : "
  cat(helptext)
  cat(paste0(defaultArgs[['SAMPLE_YML']][['description']], 
      "\n\n",
      "Optional Run-level Arguments:\n"))
  for (arg in names(defaultArgs)){
    if (arg == 'SAMPLE_YML' || !defaultArgs[[arg]][['global']]){
      next
    }
    cat(paste0(arg, 
               " : ",
               defaultArgs[[arg]][['description']],
               " [default=",
               defaultArgs[[arg]][['value']],
               "]\n"))
  }
  cat("\nOptional Sample-level Arguments:\n")
  for (arg in names(defaultArgs)){
    if (arg == 'SAMPLE_YML' || defaultArgs[[arg]][['global']]){
      next
    }
    cat(paste0(arg, 
               " : ",
               defaultArgs[[arg]][['description']],
               " [default=",
               defaultArgs[[arg]][['value']],
               "]\n"))
  }

  cat(paste0("\nThe format of SAMPLE_YML is :\n", 
            "---\n",
            "SAMPLE_NAME1:\n",
            "    GEX_datadir: /absolute/path/to/GEX_cellranger_outs/\n",
            "    cells_loaded: 25000\n",
            "    ADT_datadir: /absolute/path/to/ADT_cellranger_outs/\n",
            "    IDX_datadir: /absolute/path/to/IDX_cellranger_outs/\n",
            "    BCR_datadir: /absolute/path/to/BCR_cellranger_outs/\n",
            "    TCR_datadir: /absolute/path/to/TCR_cellranger_outs/\n",
            "    multiplet_rate: 0.35\n",
            "    metadata:\n",
            "        metadata_1_name: /absolute/path/to/tsv_of_metadata_1\n",
            "        metadata_1_name: /absolute/path/to/tsv_of_metadata_2\n",
            "    OPTIONAL_ARGUMENT_X: VALUE\n",
            "SAMPLE_NAME2:\n",
            "    GEX_datadir: /absolute/path/to/GEX_cellranger_outs/\n",
            "    cells_loaded: 15000\n",
            "...\n\n",
            "Where:\n",
            "\t1. GEX_datadir and cells_loaded are required entry but the others are optional.\n",
            "\t2. cells_loaded is expected to be an integer.\n",
            "\t3. The XXX_cellranger_outs (XXX = GEX, ADT, IDX) dirs will contain a folder named \n",
            "\t   raw_feature_bc_matrix that contains the raw counts matix.\n",
            "\t4. The XXX_cellranger_outs (XXX = TCR, BCR) dirs will contain a files named clonotypes.csv\n",
            "\t   and filtered_contig_annotations.csv.\n",
            "\t5. metadata include tsvs of metadata that need to be injected into the Seurat object.\n",
            "\t   The first column of each tsv is expected to be cell barcodes.\n",
            "\t6. multiplet_rate is a fraction between 0 and 1. If not supplied it is estimated based on\n",
            "\t   data from the 10x knowledgebase. (https://kb.10xgenomics.com/hc/en-us/articles/360001378811)\n",
            "\t7. OPTIONAL_ARGUMENT_X can be any optional sample-level argument described above and \n",
            "\t   the value in `SAMPLE_YML` will have higher precenence than the global value.\n\n",
            "NOTE:\n",
            "1. If `OUT_FOLDER` is not provided for a sample in in `SAMPLE_YML`, the deafult behaviour is\n",
            "   to use `OUT_DIR/SAMPLE_NAME`, else it is used as-is.\n",
            "2. doublet detection using DoubletFinder is not yet implemented.\n\n"))
  quit(save="no")
}

user_args = commandArgs(trailingOnly=TRUE)

args <- list()
for (i in user_args){
  key_val = strsplit(i, split="=", fixed=TRUE)[[1]]
  if (length(key_val) != 2) {
    if (length(key_val) == 1 & key_val %in% c('-h', '--help')){
        print_help(defaultArgs)
      } else {
        stop(paste0("Invalid option `", i, "`. Args must be of the form KEY=VAL."))
      }
  }
  if (!key_val[1] %in% names(defaultArgs)){
    print(paste0("WARNING: Invalid key provided : ", i, ""))
    next
  } else {
    if (key_val[1] %in% names(args)){
      args[[key_val[1]]] = c(args[[key_val[1]]], defaultArgs[[key_val[1]]][['class']](key_val[2]))
    } else {
      args[[key_val[1]]] = defaultArgs[[key_val[1]]][['class']](key_val[2])
    }
  }
}

for (i in names(defaultArgs)) {
  if (i == 'SAMPLE_YML') {
    next
  }
  if (!i %in% names(args)){
    args[[i]] <- defaultArgs[[i]][['value']]
  }
}

cat("Running with CLI arguments:\n")
cat(paste0("SAMPLE_YML                : ", args$SAMPLE_YML, '\n'))
cat(paste0("WORKING_FOLDER            : ", args$WORKING_FOLDER, '\n'))
cat(paste0("RSCRIPTS_DIR              : ", args$RSCRIPTS_DIR, '\n'))
cat(paste0("GENESET_DIR               : ", args$GENESET_DIR, '\n'))
cat(paste0("OUT_FOLDER                : ", args$OUT_FOLDER, '\n'))
cat(paste0("SPECIES                   : ", args$SPECIES, '\n'))
cat(paste0("RERUN_STAGE               : ", args$RERUN_STAGE, '\n'))
cat(paste0("SECONDARY_ASSAY_NAME      : ", args$SECONDARY_ASSAY_NAME, '\n'))
cat(paste0("IDX_ASSAY_NAME            : ", args$IDX_ASSAY_NAME, '\n'))
cat(paste0("ADT_ASSAY_NAME            : ", args$ADT_ASSAY_NAME, '\n'))
cat(paste0("IDX_NORMALIZATION_METHOD  : ", args$IDX_NORMALIZATION_METHOD, '\n'))
cat(paste0("ADT_NORMALIZATION_METHOD  : ", args$ADT_NORMALIZATION_METHOD, '\n'))
cat("\n")

if (is.null(args$SAMPLE_YML)) {
  print_help(defaultArgs)
} else if (!file.exists(args$SAMPLE_YML)){
  stop("SAMPLE_YML did not exist.")
} else {
  samples <- read_yaml(args$SAMPLE_YML)
}

species_args <- list(
  human=list(
    vars_to_regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
    mito_regex = "^MT-",
    reference_dir = paste(args$GENESET_DIR,"GRCh38", sep="/"),
    CD45 = "PTPRC"
    ),
  mouse=list(
    vars_to_regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
    mito_regex = "^Mt-",
    reference_dir = paste(args$GENESET_DIR,"GRCm38", sep="/"),
    CD45 = "Ptprc"
    )
  ) 

setwd(args$WORKING_FOLDER)

suppressPackageStartupMessages({
  source(paste(args$RSCRIPTS_DIR, "identify_hto_clusters.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "generate_profile_plot.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "demux_HTOs_by_inflexions.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "parse_10x_vdj_outs.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "TriplePlot.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "get_10x_multiplet_rate.R", sep="/"))
})

cat('Populating and validating sample run parameters.\n')
for (s in names(samples)){
  cat(paste0('Validating sample `', s, '`.\n'))
  if (!all(c('GEX_datadir', 'cells_loaded') %in% names(samples[[s]]))) {
    stop(paste0("Sample `", s, "` in `SAMPLE_YML` did not have both `GEX_datadir` and `cells_loaded`"))
  }
  for (k in names(samples[[s]])){
    if (k %in% names(allowedSampleLevelArgs)) {
      samples[[s]][[k]] <- allowedSampleLevelArgs[[k]](samples[[s]][[k]])
    } else if (k %in% local_args) {
      samples[[s]][[k]] <- defaultArgs[[k]][['class']](samples[[s]][[k]])
    } else if (k %in% global_args){
      cat(paste0('WARNING: Ignoring provided run-wide parameter `', k, '`.\n'))
      samples[[s]][[k]] <- NULL
      next
    } else {
      stop(paste0('Received unknown argument in `SAMPLE_YML`.', k))
    }
  }

  if (is.null(samples[[s]]$OUT_FOLDER)) {
    if (is.na(args$OUT_FOLDER)){
      stop("`OUT_FOLDER` must be provided in `SAMPLE_YML` if it is not provided in the command line arguments.")
    }
    suppressWarnings({
      samples[[s]]$OUT_FOLDER <- normalizePath(file.path(args$OUT_FOLDER, s))
      })
  }

  for (k in local_args){
    if (! k %in% names(samples[[s]])) {
      samples[[s]][[k]] <- args[[k]]
    }
  }

  # Check SPECIES is correct for all
  if (!samples[[s]]$SPECIES %in% c('human', 'mouse')) {
    stop("SPECIES must be one of `human` or `mouse`", call.=FALSE)
  }

  # Check RERUN_STAGE is correct for all
  if (!samples[[s]]$RERUN_STAGE %in% c(TRUE, FALSE)){
    stop("RERUN_STAGE must be one of `TRUE` or `FALSE`", call.=FALSE)
  }

  # Check cells_loaded and multiplet_rate are ok
  if (!'multiplet_rate' %in% names(samples[[s]])){
    samples[[s]]$multiplet_rate <- get_10x_multiplet_rate(samples[[s]]$cells_loaded, 
                                                          precision=3, fraction=TRUE)
  }  
}
cat('`SAMPLE_YML` looks OK!\n\n')


remove_rplots <- TRUE
if (file.exists("Rplots.pdf")){
  # Don't delete an existing file. Just delete if we create it.
  remove_rplots <- FALSE
}

sobjs <- list()
for (s in names(samples)){
  cat(paste0("Processing sample `", s, "`.\n"))
  if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed_processed.RData")) | samples[[s]]$RERUN_STAGE) {
    if (file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed_processed.RData"))) {
      cat(paste0("Found `", s, "_scTransformed_processed.RData` but regenerating it as per request.\n"))
      samples[[s]]$RERUN_STAGE <- FALSE
    }
    if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed.RData")) | samples[[s]]$RERUN_STAGE) {
      if (file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed.RData"))) {
        cat(paste0("Found `", s, "_scTransformed.RData` but regenerating it as per request.\n"))
        samples[[s]]$RERUN_STAGE <- FALSE
      }
      if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered_singlets.RData")) | samples[[s]]$RERUN_STAGE) {
        if (file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered_singlets.RData"))) {
          cat(paste0("Found `", s, "_filtered_singlets.RData` but regenerating it as per request.\n"))
          samples[[s]]$RERUN_STAGE <- FALSE
        }
        if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered.RData")) | samples[[s]]$RERUN_STAGE) {
          if (file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered.RData"))) {
            cat(paste0("Found `", s, "_filtered.RData` but regenerating it as per request.\n"))
            samples[[s]]$RERUN_STAGE <- FALSE
          }
          if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_raw.RData")) | samples[[s]]$RERUN_STAGE) {
            if (file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_raw.RData"))) {
              cat(paste0("Found `", s, "_raw.RData` but regenerating it as per request.\n"))
              samples[[s]]$RERUN_STAGE <- FALSE
            }
            cat(paste0("Generating `", s, "_raw.RData` ...\n"))
            GEX_datadir = file.path(samples[[s]][['GEX_datadir']], "raw_feature_bc_matrix")
            if (!dir.exists(GEX_datadir)) {
              cat(paste0("ERROR: Could not find a directory within `GEX_datadir` for (", s, "). ",
                         "Tried to access (", GEX_datadir, ").\n"))
              next
            }

            metadata <- list()
            metadata[['tcr']] = data.frame()
            if (!is.null(samples[[s]][['TCR_datadir']])){
              cat(paste0('LOG: Reading `TCR_datadir` ...\n'))
              tcrdir <- file.path(samples[[s]][['TCR_datadir']])
              if(dir.exists(tcrdir)){
                  cat(paste0("Reading TCR information from ", tcrdir, '\n'))
                  metadata[['tcr']] <- parse_tcr_clonotype(tcrdir)
                } else {
                  cat(paste0("`TCR_datadir` specified for ", s, " was not found at the ",
                             "specified location ", tcrdir, '\n'))
                }
            } else {
              cat(paste0('LOG: No `TCR_datadir` specified for this sample.\n'))
            }

            metadata[['bcr']] = data.frame()
            if (!is.null(samples[[s]][['BCR_datadir']])){
              cat(paste0('LOG: Reading `BCR_datadir` ...\n'))
              bcrdir <- file.path(samples[[s]][['BCR_datadir']])
              if(dir.exists(bcrdir)){
                  cat(paste0("Reading BCR information from ", bcrdir, '\n'))
                  metadata[['bcr']] <- parse_bcr_clonotype(bcrdir)
                } else {
                  cat(paste0("`BCR_datadir` specified for ", s, " was not found at the ",
                             "specified location ", bcrdir, '\n'))
                }
            } else {
              cat(paste0('LOG: No `BCR_datadir` specified for this sample.\n'))
            }

            for (md in names(samples[[s]][['metadata']])){
              fn = file.path(samples[[s]][['metadata']][[md]])
              if (file.exists(fn)){
                cat(paste0('Reading metadata file: ' ,md, '\n'))
                metadata[[md]] <- read.table(fn,
                                             sep = "\t",
                                             header = TRUE,
                                             row.names=1,
                                             stringsAsFactors = TRUE)
                rownames(metadata[[md]]) <- gsub("-1$", '', rownames(metadata[[md]]))
                if (dim(metadata[[md]])[2] < 1){
                  stop(paste0('Something wrong with metadata file ', md, '. Found 0 columns... Is the file a tsv?'))
                }
              } else {
                cat(paste0("`metadata` file ", md, " specified for ", s, 
                           " was not found at the specified location ", 
                           fn, '\n'))
              }
            }

            merged_metadata <- data.frame()
            for (md in names(metadata)){
              merged_metadata <- merge(merged_metadata, metadata[[md]], by=0, all.x=TRUE, all.y=TRUE)
              rownames(merged_metadata) <- merged_metadata$Row.names
              merged_metadata$Row.names <- NULL
            }

            if (dim(merged_metadata)[1] == 0) {
              # Empty data frame. NULL is the default for metadata in a seurat object
              merged_metadata = NULL
            }

            # Create the folder if it doesn't exist
            if (!dir.exists(samples[[s]]$OUT_FOLDER)){
              dir.create(samples[[s]]$OUT_FOLDER, showWarnings=TRUE)
            }

            cat(paste0('LOG: Reading `GEX_datadir` ...\n'))
            data <- Read10X(data.dir=GEX_datadir)

            GEX_had_secondary_assay <- FALSE
            if (typeof(data) == "list"){
              cat(paste0('WARNING: Found > 1 assays in `GEX_datadir`. Please move to using 1 assay per counts matrix.\n'))
              sobjs[[s]] <- CreateSeuratObject(counts = data$`Gene Expression`,
                                               project = s,
                                               min.cells = 3,
                                               min.features = 100)

              GEX_had_secondary_assay <- TRUE
              cat(paste0('LOG: Saving secondary assay as ', samples[[s]]$SECONDARY_ASSAY_NAME, '\n'))
              secondary_assayobj <- CreateAssayObject(counts = data$`Antibody Capture`)
              sobjs[[s]][[samples[[s]]$SECONDARY_ASSAY_NAME]] <- subset(secondary_assayobj, cells=colnames(sobjs[[s]]))
              rm(secondary_assayobj)
            } else {
              sobjs[[s]] <- CreateSeuratObject(counts = data,
                                               project = s,
                                               min.cells = 3,
                                               min.features = 100)
            }
            rm(data)
            # newer versions of Seurat
            if (!is.null(merged_metadata)){
              if(endsWith(x=colnames(sobjs[[s]])[1], "-1") && 
                  !any(grepl('-1', rownames(merged_metadata)))) {
                rownames(merged_metadata) <- paste0(rownames(merged_metadata), '-1')
              }
              sobjs[[s]] <- AddMetaData(sobjs[[s]], merged_metadata)
            }

            if (!is.null(samples[[s]][['IDX_datadir']])){
              cat(paste0('LOG: Reading `IDX_datadir` ...\n'))
              if (GEX_had_secondary_assay && samples[[s]]$SECONDARY_ASSAY_NAME == samples[[s]]$IDX_ASSAY_NAME) {
                cat(paste0("ERROR: Cannot process ", s, " because the GEX_datadir had a secondary assay ",
                           "named ", samples[[s]]$SECONDARY_ASSAY_NAME, " and an IDX_datadir with the same ",
                           "requested name was provided. Rerun with a different IDX_ASSAY_NAME to continue.\n"))
                next
              }

              IDX_datadir = file.path(samples[[s]][['IDX_datadir']], "raw_feature_bc_matrix")
              if (!dir.exists(IDX_datadir)) {
                cat(paste0("ERROR: Could not find a directory within `IDX_datadir` for (", s, "). ",
                           "Tried to access (", IDX_datadir, ").\n"))
                next
              }

              data <- Read10X(data.dir=IDX_datadir)
              if (typeof(data) == "list") {
                cat(paste0("ERROR: The counts matrix for `IDX_datadir` had > 1 assays.\n"))
                next
              }
              IDX_assayobj <- CreateAssayObject(counts = data)
              sobjs[[s]][[samples[[s]]$IDX_ASSAY_NAME]] <- subset(IDX_assayobj, cells=colnames(sobjs[[s]]))
              rm(data, IDX_assayobj)
            } else {
              cat('LOG: No `IDX_datadir` specified for this sample.\n')
            }

            if (!is.null(samples[[s]][['ADT_datadir']])){
              cat(paste0('LOG: Reading `ADT_datadir` ...\n'))
              if (GEX_had_secondary_assay && samples[[s]]$SECONDARY_ASSAY_NAME == samples[[s]]$ADT_ASSAY_NAME) {
                cat(paste0("ERROR: Cannot process ", s, " because the GEX_datadir had a secondary assay ",
                           "named ", samples[[s]]$SECONDARY_ASSAY_NAME, " and an ADT_datadir with the same ",
                           "requested name was provided. Rerun with a different ADT_ASSAY_NAME to continue.\n"))
                next
              }

              ADT_datadir = file.path(samples[[s]][['ADT_datadir']], "raw_feature_bc_matrix")
              if (!dir.exists(ADT_datadir)) {
                cat(paste0("ERROR: Could not find a directory within `ADT_datadir` for (", s, "). ",
                           "Tried to access (", ADT_datadir, ").\n"))
                next
              }

              data <- Read10X(data.dir=ADT_datadir)
              if (typeof(data) == "list") {
                cat(paste0("ERROR: The counts matrix for `ADT_datadir` had > 1 assays.\n"))
                next
              }
              ADT_assayobj <- CreateAssayObject(counts = data)
              sobjs[[s]][[samples[[s]]$ADT_ASSAY_NAME]] <- subset(ADT_assayobj, cells=colnames(sobjs[[s]]))
              rm(data, ADT_assayobj)
            } else {
              cat('LOG: No `ADT_datadir` specified for this sample.\n')
            }


            # store mitochondrial percentage in object metadata
            sobjs[[s]] <- PercentageFeatureSet(sobjs[[s]],
                                               pattern = species_args[[samples[[s]]$SPECIES]]$mito_regex,
                                               col.name = "percent.mt")
            #Store ribosomal percentage in the object metadata
            ribo_genes <- read.table(paste(species_args[[samples[[s]]$SPECIES]]$reference_dir,
                                         "ribo_genes.tsv", sep="/"),
                                     sep = "\t", header=TRUE,
                                     stringsAsFactors = FALSE)
            ribo_genes <- ribo_genes[ribo_genes[["HUGO"]] %in% rownames(sobjs[[s]]), ]
            sobjs[[s]] <- PercentageFeatureSet(sobjs[[s]],
                                               features = ribo_genes[["HUGO"]],
                                               col.name = "percent.ribo")

            #Store cell cycle state in the object metadata
            cc_genes <- read.table(paste(species_args[[samples[[s]]$SPECIES]]$reference_dir,
                                         "cell_cycle_genes.tsv", sep="/"),
                                   sep = "\t", header=TRUE,
                                   stringsAsFactors = FALSE)

            if (samples[[s]]$SPECIES == "human"){
              sobjs[[s]] <- CellCycleScoring(sobjs[[s]],
                                             s.features = cc_genes[cc_genes$stage=="G1-S", "HUGO"],
                                             g2m.features = cc_genes[cc_genes$stage=="G2-M", "HUGO"],
                                             nbin = 12)
            } else if (samples[[s]]$SPECIES == "mouse") {
              sobjs[[s]] <- CellCycleScoring(sobjs[[s]],
                                             s.features = cc_genes[cc_genes$stage=="G1-S", "MGI"],
                                             g2m.features = cc_genes[cc_genes$stage=="G2-M", "MGI"],
                                             nbin = 12)
            }

            plot_all_profiles(sobjs[[s]],
                              out_prefix=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_dualscatter_pre"),
                              plot_extra=FALSE,
                              other_assay_names=unique(c(samples[[s]]$SECONDARY_ASSAY_NAME, 
                                                         samples[[s]]$IDX_ASSAY_NAME, 
                                                         samples[[s]]$ADT_ASSAY_NAME)))

            cutoffs = list()
            for (hl in c('high', 'low')){
              for (feature in c('percent.mt', 'percent.ribo')){
                cutoffs[[paste0(feature, '.', hl)]] = NA
              }
              for (assay in Assays(sobjs[[s]])){
                cutoffs[[paste0('nCount_', assay, '.', hl)]] = NA
                cutoffs[[paste0('nFeature_', assay, '.', hl)]] = NA
              }
            }

            if ("DROPLET.TYPE" %in% colnames(sobjs[[s]]@meta.data)) {
              # We have passed in Demuxlet/Freemuxlet info
              plot_all_profiles(sobjs[[s]],
                                out_prefix=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_droplet_dualscatter_pre"),
                                plot_extra=FALSE,
                                scatter_color="DROPLET.TYPE",
                                other_assay_names=unique(c(samples[[s]]$SECONDARY_ASSAY_NAME, 
                                                           samples[[s]]$IDX_ASSAY_NAME, 
                                                           samples[[s]]$ADT_ASSAY_NAME)))
              cutoffs[["DROPLET.TYPE.keep"]] = NA
              cutoffs[["BEST.GUESS.keep"]] = NA
            }

            writeLines(colnames(sobjs[[s]]), con=paste0(samples[[s]]$OUT_FOLDER, '/', "barcodes_of_interest.list"))

            write_yaml(cutoffs, file=paste0(samples[[s]]$OUT_FOLDER, '/', "cutoffs.yml"))

            assign(s, sobjs[[s]])
            save(list=s, file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_raw.RData"))
            next
          } else {
            print(paste0("Using existing `", s, "_raw.RData`."))
            load(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_raw.RData"))
            sobjs[[s]] <- get(s)
          }  # END STAGE 1
          rm(list=s)
          print(paste0("Generating `", s, "_filtered.RData` ..."))
          if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', "cutoffs.yml"))) {
            print("Cannot continue without a cutoffs.yml file in the output directory")
            next
          }
          cutoffs <- read_yaml(paste0(samples[[s]]$OUT_FOLDER, '/', "cutoffs.yml"))
          # Filter cells with high mito content (dying/dead cells)
          keep_rownames = rep(TRUE, dim(sobjs[[s]]@meta.data)[1])
          total_cells = c(dim(sobjs[[s]]@meta.data)[1], 100)
          additional_text = ""

          for (filter_key in names(cutoffs)){
            filter_feature = strsplit(filter_key, split = ".",
                                      fixed = TRUE)[[1]]
            suppmsg <- assert_that(length(filter_feature) >= 2,
                                   msg=paste0("Names in cutoffs.yml must ",
                                              "be of the form `feature.high`, ",
                                              "`feature.low`, or `feature.keep`. ",
                                              "Got ", filter_key))
            filter_cat <- filter_feature[length(filter_feature)]
            filter_feature <- paste(filter_feature[1:(length(filter_feature)-1)],
                                    collapse=".")
            suppmsg <- assert_that(filter_cat %in% c("high", "low", "keep"),
                                   msg=paste0("Names in cutoffs.yml must ",
                                              "be of the form `feature.high`, ",
                                              "`feature.low`, or `feature.keep`. ",
                                              "Got ", filter_key))
            if (!filter_feature %in% colnames(sobjs[[s]]@meta.data)){
              print(paste0("Could not find ", filter_feature, " in the metadata ",
                           "for ", s))
              next
            }
            if(!is.na(cutoffs[[filter_key]])){
              if (filter_cat == "high") {
                keep_rownames <- keep_rownames & sobjs[[s]]@meta.data[[filter_feature]] <= cutoffs[[filter_key]]
                filter_cat_text <- " > "
              } else if (filter_cat == "low") {
                keep_rownames <- keep_rownames & sobjs[[s]]@meta.data[[filter_feature]] >= cutoffs[[filter_key]]
                filter_cat_text <- " < "
              } else {
                keep_rownames <- keep_rownames &
                  sobjs[[s]]@meta.data[[filter_feature]] %in% strsplit(cutoffs[[filter_key]], ',')[[1]]
                filter_cat_text <- " not in "
              }

              print(paste0("Dropping ", additional_text,
                           total_cells[1] - sum(keep_rownames),
                           " cells (",
                           round((total_cells[1] - sum(keep_rownames))/length(keep_rownames) * 100, 2),
                           "%) for having `", filter_feature,
                           filter_cat_text,
                           cutoffs[[filter_key]], "`."))
              additional_text = "an additional "
              total_cells = c(sum(keep_rownames),
                              round(sum(keep_rownames)/length(keep_rownames) * 100, 2))
            }
          }

          print(paste0("Dropping a total of ", sum(keep_rownames==FALSE),
                       " cells (",
                       round(sum(keep_rownames==FALSE)/length(keep_rownames) *100, 2),
                       "%)."))
          print(paste0("Retaining a total of ", sum(keep_rownames),
                       " cells (",
                       round(sum(keep_rownames)/length(keep_rownames) *100, 2),
                       "%)."))

          sobjs[[s]] <- subset(sobjs[[s]],
                               cells = rownames(sobjs[[s]]@meta.data[keep_rownames, ]))
          keep_features <- c()
          for (assay in names(sobjs[[s]]@assays)){
            if (assay == "RNA"){
              next
            }
            keep_features <- c(keep_features, rownames(sobjs[[s]]@assays[[assay]]@counts))
          }
          keep_genes <- rownames(sobjs[[s]]@assays$RNA@counts)[Matrix::rowSums(sobjs[[s]]@assays$RNA@counts>0)>3]
          print(paste0("Dropping ",
                       (length(rownames(sobjs[[s]]@assays$RNA@counts)) -
                        length(keep_genes)),
                       " genes for being in fewer than 3 cells"))
          sobjs[[s]] <- subset(sobjs[[s]], features = c(keep_genes, keep_features))

          plot_all_profiles(sobjs[[s]],
                            out_prefix=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_dualscatter_post"),
                            plot_extra=TRUE,
                            other_assay_names=unique(c(samples[[s]]$SECONDARY_ASSAY_NAME, 
                                                       samples[[s]]$IDX_ASSAY_NAME, 
                                                       samples[[s]]$ADT_ASSAY_NAME)))

          if ("DROPLET.TYPE" %in% colnames(sobjs[[s]]@meta.data)) {
            # We have passed in Demuxlet/Freemuxlet info
            plot_all_profiles(sobjs[[s]],
                              out_prefix=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_droplet_dualscatter_post"),
                              plot_extra=TRUE,
                              scatter_color="DROPLET.TYPE",
                              other_assay_names=unique(c(samples[[s]]$SECONDARY_ASSAY_NAME, 
                                                         samples[[s]]$IDX_ASSAY_NAME, 
                                                         samples[[s]]$ADT_ASSAY_NAME)))
            sobjs[[s]]@meta.data$SAMPLE.by.SNPs = as.factor(sapply(as.vector(sobjs[[s]]@meta.data$BEST.GUESS), function(x) {
              temp <- paste(sort(unique(strsplit(x, ',')[[1]])), collapse='_')
              }))
          }

          # Now that we"ve removed bogus cells, let's normalize the counts (if present)
          if (samples[[s]]$IDX_ASSAY_NAME %in% names(sobjs[[s]]@assays)){
            sobjs[[s]] <- NormalizeData(sobjs[[s]], assay = samples[[s]]$IDX_ASSAY_NAME,
                                        normalization.method = samples[[s]]$IDX_NORMALIZATION_METHOD)
            sobjs[[s]] <- ScaleData(sobjs[[s]], assay = samples[[s]]$IDX_ASSAY_NAME)
          }
          if (samples[[s]]$ADT_ASSAY_NAME %in% names(sobjs[[s]]@assays)){
            sobjs[[s]] <- NormalizeData(sobjs[[s]], assay = samples[[s]]$ADT_ASSAY_NAME,
                                        normalization.method = samples[[s]]$ADT_NORMALIZATION_METHOD)
            sobjs[[s]] <- ScaleData(sobjs[[s]], assay = samples[[s]]$ADT_ASSAY_NAME)
          }
          assign(s, sobjs[[s]])
          save(list=s, file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered.RData"))
          next
        } else {
          print(paste0("Using existing `", s, "_filtered.RData`."))
          load(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered.RData"))
          sobjs[[s]] <- get(s)
        }  # END STAGE 2
        rm(list=s)
        if (samples[[s]]$IDX_ASSAY_NAME %in% names(sobjs[[s]]@assays)){
          print(paste0("Generating `", s, "_filtered_singlets.RData` ..."))
          if (!file.exists(paste0(samples[[s]]$OUT_FOLDER, '/', "IDX_map.tsv"))){
            print("Cannot continue without an IDX_map.tsv file")
            next
          }
          IDX_map <- read.table(paste0(samples[[s]]$OUT_FOLDER, '/', "IDX_map.tsv"), sep="\t",
                                header=TRUE, row.names=1, stringsAsFactors=FALSE)
          suppmgg <- assert_that(all('sample_name' %in% colnames(IDX_map)),
                                 msg="IDX_map.tsv must have a column named `sample_name`")
          # Hacky but it works ¯\_(ツ)_/¯
          IDX_map <- sapply(colnames(IDX_map), function(x) {
                              y = IDX_map[,x, drop=T]
                              names(y) = rownames(IDX_map)
                              y
                              }, simplify = FALSE, USE.NAMES = TRUE)
          sample_names <- names(IDX_map$sample_name)
          sample_names <- sample_names[!sample_names %in% c("NEGATIVE", "MULTIPLET", "OTHER")]
          suppmgg <- assert_that(all(sample_names %in% rownames(sobjs[[s]]@assays[[samples[[s]]$IDX_ASSAY_NAME]])),
                                 msg="Not all hashtags in IDX_map.tsv are in the seurat object")
          keep_features <- c(rownames(sobjs[[s]]@assays$RNA),
                             sample_names)
          sobjs[[s]] <- subset(sobjs[[s]], features = keep_features)

          if ("threshold" %in% names(IDX_map)){
            inflexions <- IDX_map$threshold[sample_names]
            if (any(is.na(inflexions))){
              inflexions <- NULL
            } else {
              print("Using user-provided thresholds for demultiplexing.")
            }
          } else {
            inflexions <- NULL
          }

          if ("color" %in% names(IDX_map)){
            sample_colors <- IDX_map$color
            suppmsg <- assert_that(all(c("NEGATIVE", "MULTIPLET", "OTHER") %in% names(sample_colors)),
                                   msg=paste0("Cannot use user-provided colors if values ",
                                              "for NEGATIVE,  MULTIPLET, and OTHER are not ",
                                              "provided"))
            if (any(is.na(sample_colors))){
              sample_colors <- NULL
            } else {
              #TODO
              sample_colors <- NULL
              print("Can't use user-provided colors for plotting yet.")
              #print("Using user-provided colors for plotting.")
            }
          } else {
            sample_colors <- NULL
          }

          demux_results <- demux_by_inflexions(sobjs[[s]],
                                               sample_names=IDX_map$sample_name[sample_names],
                                               inflexions=inflexions,
                                               sample_colors=sample_colors,
                                               assay_name=samples[[s]]$IDX_ASSAY_NAME)

          nplots <- length(demux_results$background_plots)

          png(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_raw_hto_pairplots.png"), width=nplots*750,
              height=nplots*750, units="px")
          print(plot_grid(plotlist=unlist(demux_results[['raw_plots']], recursive = FALSE),
                          ncol=nplots))
          dev.off()
          png(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scaled_hto_pairplots.png"), width=nplots*750,
              height=nplots*750, units="px")
          print(plot_grid(plotlist=unlist(demux_results[['scaled_plots']], recursive=FALSE),
                          ncol=nplots))
          dev.off()
          png(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_hto_background_profiles.png"), width=1500,
              height=750*max(1, floor(nplots/2)), units="px")
          print(plot_grid(plotlist=demux_results[['background_plots']], ncol=2))
          dev.off()

          sobjs[[s]]@meta.data$SAMPLE.by.ABs <- demux_results$cell_ids
          sobjs[[s]]@meta.data$DROPLET.TYPE.by.ABs <- demux_results$droplet_type

          IDX_metadata <- sobjs[[s]]@meta.data[, c('SAMPLE.by.ABs', 'DROPLET.TYPE.by.ABs')]
          write.table(IDX_metadata,
                      file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered_IDX_calls.tsv"),
                      quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

          drf <- data.frame(sample_name=unlist(IDX_map$sample_name[sample_names]),
                  threshold=demux_results[['inflexions']])
          write.table(drf, file=paste0(samples[[s]]$OUT_FOLDER, '/', "IDX_map_generated.tsv"),
                      quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

          capture.output(demux_results[['secondary_inflexions']],
                         file=paste0(samples[[s]]$OUT_FOLDER, '/', "IDX_additional_inflexions.txt"))

          write.table(demux_results[["summary"]], file=paste0(samples[[s]]$OUT_FOLDER, '/', "IDX_summary.tsv"),
                      quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

          sobjs[[s]]$temp <- sobjs[[s]]$SAMPLE.by.ABs
          sobjs[[s]]$temp[sobjs[[s]]$DROPLET.TYPE.by.ABs == 'MULTIPLET'] = 'MULTIPLET'
          sobjs[[s]]$temp <- factor(sobjs[[s]]$temp, levels = c('NEGATIVE', 
                                                                unname(sort(IDX_map$sample_name)),
                                                                'MULTIPLET'))
          Idents(sobjs[[s]]) <- sobjs[[s]]$temp
          sobjs[[s]]$temp <- NULL
          png(paste0(samples[[s]]$OUT_FOLDER, '/', s, '_ridgeplot.png'), width=1500,
              height=ceiling(length(sample_names)/2)*750, units = 'px')
          print(RidgePlot(sobjs[[s]],
                          assay=samples[[s]]$IDX_ASSAY_NAME,
                          features=sample_names,
                          ncol = 2))
          dev.off()
          sobjs[[s]] <- subset(sobjs[[s]],
                               cells=colnames(sobjs[[s]])[sobjs[[s]]$DROPLET.TYPE.by.ABs == 'SINGLET'])
          # This will just become SINGLET. drop it to reduce object size
          sobjs[[s]]$DROPLET.TYPE.by.ABs <- NULL 
          sobjs[[s]]$SAMPLE.by.ABs <- factor(sobjs[[s]]$SAMPLE.by.ABs,
                                             levels=unname(sort(IDX_map$sample_name)))
          assign(s, sobjs[[s]])
          save(list=s, file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered_singlets.RData"))
          next
        } else {
          print(paste0("Could not find an assay named ", samples[[s]]$IDX_ASSAY_NAME,
                       " in the object. Continuing on to run SCTransform. ",
                       "If you have multiplexed this data and have stored the ",
                       "counts in another assay name, delete the SCTransformed ",
                       "RData rerun this stage with IDX_ASSAY_NAME=XXX."))
        }
      } else {
        print(paste0("Using existing `", s, "_filtered_singlets.RData`."))
        load(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_filtered_singlets.RData"))
        sobjs[[s]] <- get(s)
      }  # END STAGE 3
      rm(list=s)
      print(paste0("Generating `", s, "_scTransformed.RData` ..."))
      sobjs[[s]] <- SCTransform(sobjs[[s]],
                                vars.to.regress = species_args[[samples[[s]]$SPECIES]]$vars_to_regress,
                                return.only.var.genes = FALSE,
                                verbose = FALSE)

      # Get the Principal components for the object
      sobjs[[s]] <- RunPCA(sobjs[[s]], verbose = FALSE)
      # Verify the confounding vars have been regressed out
      pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_confound_pca.pdf"))
      print(FeaturePlot(sobjs[[s]], reduction = "pca",
                        features=species_args[[samples[[s]]$SPECIES]]$vars_to_regress))
      dev.off()
      assign(s, sobjs[[s]])
      save(list=s, file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed.RData"))
      next
    } else {
      print(paste0("Using existing `", s, "_scTransformed.RData`."))
      load(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed.RData"))
      sobjs[[s]] <- get(s)
    }  # END STAGE 4
    rm(list=s)
    print(paste0("Generating `", s, "_scTransformed_processed.RData`..."))
    sobjs[[s]] <- RunUMAP(sobjs[[s]],
                     dims = 1:30,  # Num PCs to use
                     n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                     min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                     spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                     a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                     b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                     verbose = FALSE)
    pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_umap.pdf"))
    print(DimPlot(sobjs[[s]]))
    dev.off()

    # Verify the confounding vars have been regressed out
    pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_confound_umap.pdf"))
    print(FeaturePlot(sobjs[[s]], features=species_args[[samples[[s]]$SPECIES]]$vars_to_regress))
    dev.off()
    # Calculate the neighborhood graph
    sobjs[[s]] <- FindNeighbors(sobjs[[s]],
                          dims = 1:30,  # Num PCs to use
                          k.param = 20,  # k for the knn algorithm
                          verbose = FALSE)
    # Use the neighborhood graph to cluster the data
    sobjs[[s]] <- FindClusters(sobjs[[s]], verbose = FALSE,
                               algorithm = 4)  # Use Leiden

    # View the cluters UMAP
    pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_umap.pdf"))
    print(DimPlot(sobjs[[s]], label=TRUE))
    dev.off()

    if ("SAMPLE.by.SNPs" %in% colnames(sobjs[[s]]@meta.data)) {
      # We have passed in Demuxlet/Freemuxlet info
      pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_samples_by_SNPs_umap.pdf"))
      print(DimPlot(sobjs[[s]], group.by="SAMPLE.by.SNPs"))
      dev.off()
    }

    if ("SAMPLE.by.ABs" %in% colnames(sobjs[[s]]@meta.data)){
      # View the clusters UMAP
      pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_samples_by_ABs_umap.pdf"))
      print(DimPlot(sobjs[[s]], group.by="SAMPLE.by.ABs"))
      dev.off()
    }

    pdf(paste0(samples[[s]]$OUT_FOLDER, '/', s, "_immune_triplot_umap.pdf"))
    tryCatch(expr=print(TriPlot(sobjs[[s]], features=species_args[[samples[[s]]$SPECIES]]$CD45, reduction.use="umap", group.by="seurat_clusters")),
             error=function(e) {print(paste0("WARNING: Could not print CD45(",
                                species_args[[samples[[s]]$SPECIES]]$CD45,
                                ") expression since it was not found in the object"))
                                ggplot() + theme_void()},
             finally=function(e){dev.off()})

    # TODO: Add doubletfinder code

    assign(s, sobjs[[s]])
    save(list=s, file=paste0(samples[[s]]$OUT_FOLDER, '/', s, "_scTransformed_processed.RData"))
    next
  } else {
    print(paste0("Nothing more to do for ", s, " since `", s,
                 "_scTransformed_processed.RData` exists."))
  }
}  # END STAGE 5

if (remove_rplots & file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

uuid <- paste(sample(c(LETTERS,0:9), 10, replace = T), collapse="")
write(paste(names(warnings()), collapse="\n"), file=paste0("run_warnings_", uuid, ".txt"))
cat(paste0("Wrote run warnings to run_warnings_", uuid, ".txt\n"))
