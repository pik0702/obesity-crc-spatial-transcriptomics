#!/usr/bin/env Rscript
Sys.setenv(
  GEOMX_DATA_DIR = "C:/Users/user/Documents/kim"
)
###############################################################################
# GeoMx DSP preprocessing, QC, gene filtering, and Q3 normalization
#
# Study:
#   Obesity-associated colorectal cancer spatial transcriptomics
#
# Purpose:
#   1. Import GeoMx DCC/PKC files and de-identified ROI annotation.
#   2. Perform segment-level and probe-level QC.
#   3. Aggregate probe counts to gene/target level.
#   4. Calculate LOQ-based gene detection rates.
#   5. Retain genes detected above LOQ in >=5% of ROIs
#      (negative control targets are retained internally for normalization/QC).
#   6. Perform Q3 (75th percentile) normalization.
#   7. Export filtered gene lists, ROI metadata, raw counts,
#      Q3-normalized expression, and log2(Q3 + 1) expression.
#
# Input directory structure:
#
#   data/GeoMx/
#   ├── dccs/
#   │   └── *.dcc
#   ├── pkcs/
#   │   └── *.pkc
#   └── annotation/
#       └── *.xlsx
#
# Usage in RStudio:
#   Edit DATA_DIR if needed, then source this script.
#
# Usage from command line:
#   Rscript GeoMx/01_preprocessing_QC_normalization.R
#
###############################################################################


###############################################################################
# 0. Configuration
###############################################################################

DATA_DIR <- Sys.getenv(
  "GEOMX_DATA_DIR",
  unset = file.path("data", "GeoMx")
)

OUTPUT_DIR <- Sys.getenv(
  "GEOMX_OUTPUT_DIR",
  unset = file.path(
    "results",
    "GeoMx",
    "01_preprocessing_QC_normalization"
  )
)

ANNOTATION_SHEET <- "Template"
ANNOTATION_DCC_COLUMN <- "Sample_ID"

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
###############################################################################
# 1. Required packages
###############################################################################

required_packages <- c(
  "NanoStringNCTools",
  "GeomxTools",
  "GeoMxWorkflows",
  "Biobase"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall the required packages before running this script."
  )
}

suppressPackageStartupMessages({
  library(NanoStringNCTools)
  library(GeomxTools)
  library(GeoMxWorkflows)
  library(Biobase)
})


###############################################################################
# 2. Locate input files
###############################################################################

DCCFiles <- list.files(
  file.path(DATA_DIR, "dccs"),
  pattern = "\\.dcc$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

PKCFiles <- list.files(
  file.path(DATA_DIR, "pkcs"),
  pattern = "\\.pkc$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

SampleAnnotationFile <- list.files(
  file.path(DATA_DIR, "annotation"),
  pattern = "\\.xlsx$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# Exclude temporary Excel lock files
SampleAnnotationFile <- SampleAnnotationFile[
  !grepl("^~\\$", basename(SampleAnnotationFile))
]



if (length(DCCFiles) == 0) {
  stop("No DCC files were found under: ", file.path(DATA_DIR, "dccs"))
}

if (length(PKCFiles) == 0) {
  stop("No PKC files were found under: ", file.path(DATA_DIR, "pkcs"))
}

if (length(SampleAnnotationFile) != 1) {
  stop(
    "Exactly one non-temporary annotation .xlsx file is required under ",
    file.path(DATA_DIR, "annotation"),
    ". Found: ",
    length(SampleAnnotationFile)
  )
}

message("DCC files: ", length(DCCFiles))
message("PKC files: ", length(PKCFiles))
message("Annotation: ", SampleAnnotationFile)


###############################################################################
# 3. Import GeoMx data
###############################################################################

geomx <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = ANNOTATION_SHEET,
  phenoDataDccColName = ANNOTATION_DCC_COLUMN,
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

# Shift zero counts to one according to the GeoMx workflow.
geomx <- shiftCountsOne(
  geomx,
  useDALogic = TRUE
)

pkcs <- annotation(geomx)

modules <- unique(
  gsub(
    "\\.pkc$",
    "",
    basename(pkcs),
    ignore.case = TRUE
  )
)

message(
  "Imported object: ",
  nrow(geomx), " probes x ", ncol(geomx), " ROIs"
)


###############################################################################
# 4. Segment-level QC
###############################################################################

# IMPORTANT:
# These numerical cutoffs are preserved from the analysis script used for
# this study. Verify them against the final Methods before public release.
QC_params <- list(
  minSegmentReads = 1000,
  percentTrimmed = 80,
  percentStitched = 50,
  percentAligned = 50,
  percentSaturation = 50,
  minNegativeCount = 1,
  maxNTCCount = 100000,
  minNuclei = 1,
  minArea = 100
)

geomx <- setSegmentQCFlags(
  geomx,
  qcCutoffs = QC_params
)

QCResults <- protocolData(geomx)[["QCFlags"]]
flag_columns <- colnames(QCResults)

QC_Summary <- data.frame(
  QC_metric = flag_columns,
  Pass = colSums(!QCResults[, flag_columns, drop = FALSE]),
  Warning = colSums(QCResults[, flag_columns, drop = FALSE]),
  row.names = NULL,
  check.names = FALSE
)

QCResults$QCStatus <- apply(
  QCResults[, flag_columns, drop = FALSE],
  1L,
  function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  }
)

QC_total <- data.frame(
  QC_metric = "TOTAL_ROIS",
  Pass = sum(QCResults$QCStatus == "PASS"),
  Warning = sum(QCResults$QCStatus == "WARNING"),
  row.names = NULL
)

QC_Summary <- rbind(QC_Summary, QC_total)

write.csv(
  QC_Summary,
  file.path(OUTPUT_DIR, "GeoMx_segment_QC_summary.csv"),
  row.names = FALSE
)

# Save per-ROI QC flags before removing failed ROIs.
QC_per_ROI <- cbind(
  ROI = rownames(QCResults),
  as.data.frame(QCResults, check.names = FALSE)
)

write.csv(
  QC_per_ROI,
  file.path(OUTPUT_DIR, "GeoMx_segment_QC_flags.csv"),
  row.names = FALSE
)


###############################################################################
# 5. Negative-control statistics and removal of failed ROIs
###############################################################################

negativeGeoMeans <- esBy(
  negativeControlSubset(geomx),
  GROUP = "Module",
  FUN = function(x) {
    assayDataApply(
      x,
      MARGIN = 2,
      FUN = ngeoMean,
      elt = "exprs"
    )
  }
)

protocolData(geomx)[["NegGeoMean"]] <- negativeGeoMeans

negCols <- paste0(
  "NegGeoMean_",
  modules
)

# Preserve the workflow behavior used in the original analysis.
if (length(negCols) > 0) {
  neg_means_from_sData <- sData(geomx)[["NegGeoMean"]]
  
  if (!is.null(neg_means_from_sData)) {
    pData(geomx)[, negCols] <- neg_means_from_sData
  }
}

# Remove failed ROIs.
geomx <- geomx[
  ,
  QCResults$QCStatus == "PASS"
]

message(
  "ROIs retained after segment QC: ",
  ncol(geomx)
)


###############################################################################
# 6. Probe-level QC
###############################################################################

geomx <- setBioProbeQCFlags(
  geomx,
  qcCutoffs = list(
    minProbeRatio = 0.1,
    percentFailGrubbs = 20
  ),
  removeLocalOutliers = TRUE
)

ProbeQCResults <- fData(geomx)[["QCFlags"]]

ProbeQC_summary <- data.frame(
  Passed = sum(rowSums(ProbeQCResults[, -1, drop = FALSE]) == 0),
  Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
  Local = sum(
    rowSums(ProbeQCResults[, -c(1, 2), drop = FALSE]) > 0 &
      !ProbeQCResults$GlobalGrubbsOutlier
  )
)

write.csv(
  ProbeQC_summary,
  file.path(OUTPUT_DIR, "GeoMx_probe_QC_summary.csv"),
  row.names = FALSE
)

# Remove probes failing LowProbeRatio or GlobalGrubbsOutlier QC.
ProbeQCPassed <- subset(
  geomx,
  fData(geomx)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
    fData(geomx)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE
)

geomx <- ProbeQCPassed

message(
  "Probes retained after probe QC: ",
  nrow(geomx)
)


###############################################################################
# 7. Aggregate probe counts to gene/target level
###############################################################################

target_geomx <- aggregateCounts(geomx)

message(
  "After target aggregation: ",
  nrow(target_geomx), " targets x ",
  ncol(target_geomx), " ROIs"
)


###############################################################################
# 8. LOQ calculation and gene detection rate
###############################################################################

LOQ_SD_CUTOFF <- 2
MIN_LOQ <- 2

LOQ <- data.frame(
  row.names = colnames(target_geomx)
)

usable_modules <- character(0)

for (module in modules) {
  
  vars <- paste0(
    c("NegGeoMean_", "NegGeoSD_"),
    module
  )
  
  if (all(vars %in% colnames(pData(target_geomx)))) {
    
    LOQ[, module] <- pmax(
      MIN_LOQ,
      pData(target_geomx)[, vars[1]] *
        pData(target_geomx)[, vars[2]] ^ LOQ_SD_CUTOFF
    )
    
    usable_modules <- c(
      usable_modules,
      module
    )
  }
}

usable_modules <- unique(usable_modules)

if (length(usable_modules) == 0) {
  stop(
    paste0(
      "LOQ could not be calculated because the expected NegGeoMean_/NegGeoSD_ ",
      "columns were not found in the target-level object. Inspect colnames(pData(target_geomx))."
    )
  )
}

pData(target_geomx)$LOQ <- LOQ

LOQ_Mat <- NULL

for (module in usable_modules) {
  
  ind <- fData(target_geomx)$Module == module
  
  if (!any(ind)) {
    next
  }
  
  Mat_i <- t(
    esApply(
      target_geomx[ind, ],
      MARGIN = 1,
      FUN = function(x) {
        x > LOQ[, module]
      }
    )
  )
  
  LOQ_Mat <- rbind(
    LOQ_Mat,
    Mat_i
  )
}

target_names <- as.character(
  fData(target_geomx)$TargetName
)

match_idx <- match(
  target_names,
  rownames(LOQ_Mat)
)

if (anyNA(match_idx)) {
  stop(
    "Some aggregated targets could not be matched to the LOQ detection matrix."
  )
}

LOQ_Mat <- LOQ_Mat[
  match_idx,
  ,
  drop = FALSE
]

rownames(LOQ_Mat) <- target_names

fData(target_geomx)$DetectedSegments <- rowSums(
  LOQ_Mat,
  na.rm = TRUE
)

fData(target_geomx)$DetectionRate <- (
  fData(target_geomx)$DetectedSegments /
    ncol(target_geomx)
)


###############################################################################
# 9. Gene filtering
###############################################################################

GENE_DETECTION_THRESHOLD <- 0.05

negative_idx <- (
  as.character(fData(target_geomx)$CodeClass) == "Negative"
)

keep_gene <- (
  fData(target_geomx)$DetectionRate >= GENE_DETECTION_THRESHOLD |
    negative_idx
)

target_geomx <- target_geomx[
  keep_gene,
  ,
  drop = FALSE
]

message(
  "Targets retained after >= ",
  GENE_DETECTION_THRESHOLD * 100,
  "% ROI detection filter: ",
  nrow(target_geomx)
)

# Biological genes only (negative controls removed for exported expression data).
biological_idx <- (
  as.character(fData(target_geomx)$CodeClass) != "Negative"
)

biological_target_names <- as.character(
  fData(target_geomx)$TargetName[biological_idx]
)

message(
  "Biological genes retained: ",
  length(biological_target_names)
)


###############################################################################
# 10. Q3 normalization
###############################################################################

target_geomx <- normalize(
  target_geomx,
  norm_method = "quant",
  desiredQuantile = 0.75,
  toElt = "q_norm"
)

raw_count_matrix <- assayDataElement(
  target_geomx,
  elt = "exprs"
)

q3_matrix <- assayDataElement(
  target_geomx,
  elt = "q_norm"
)

# Export only biological genes for the processed expression matrices.
raw_count_bio <- raw_count_matrix[
  biological_idx,
  ,
  drop = FALSE
]

q3_bio <- q3_matrix[
  biological_idx,
  ,
  drop = FALSE
]

rownames(raw_count_bio) <- biological_target_names
rownames(q3_bio) <- biological_target_names

log2_q3_bio <- log2(
  q3_bio + 1
)


###############################################################################
# 11. Export results
###############################################################################

# Helper: write a gene-by-ROI expression matrix with a Gene column.
write_expression_csv <- function(mat, filename) {
  
  out <- data.frame(
    Gene = rownames(mat),
    as.data.frame(
      mat,
      check.names = FALSE
    ),
    check.names = FALSE
  )
  
  write.csv(
    out,
    file.path(OUTPUT_DIR, filename),
    row.names = FALSE,
    quote = FALSE
  )
}


# 11.1 Filtered gene list + detection metadata
gene_metadata <- as.data.frame(
  fData(target_geomx),
  check.names = FALSE
)

gene_metadata$TargetName <- as.character(
  gene_metadata$TargetName
)

gene_metadata_bio <- gene_metadata[
  as.character(gene_metadata$CodeClass) != "Negative",
  ,
  drop = FALSE
]

write.csv(
  gene_metadata_bio,
  file.path(
    OUTPUT_DIR,
    "GeoMx_filtered_gene_metadata.csv"
  ),
  row.names = FALSE
)

writeLines(
  biological_target_names,
  con = file.path(
    OUTPUT_DIR,
    "GeoMx_filtered_gene_list.txt"
  )
)


# 11.2 Raw expression matrix after QC/gene filtering
write_expression_csv(
  raw_count_bio,
  "GeoMx_filtered_raw_expression.csv"
)


# 11.3 Q3-normalized expression matrix
write_expression_csv(
  q3_bio,
  "GeoMx_Q3_normalized_expression.csv"
)


# 11.4 log2(Q3 + 1) expression matrix
write_expression_csv(
  log2_q3_bio,
  "GeoMx_log2_Q3_normalized_expression.csv"
)


# 11.5 ROI metadata after QC
roi_metadata <- data.frame(
  ROI = rownames(pData(target_geomx)),
  as.data.frame(
    pData(target_geomx),
    check.names = FALSE
  ),
  check.names = FALSE
)

write.csv(
  roi_metadata,
  file.path(
    OUTPUT_DIR,
    "GeoMx_ROI_metadata_after_QC.csv"
  ),
  row.names = FALSE
)


# 11.6 Save the processed GeoMx object for downstream scripts.
saveRDS(
  target_geomx,
  file.path(
    OUTPUT_DIR,
    "GeoMx_preprocessed_Q3_normalized.rds"
  )
)


# 11.7 Record software versions.
writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)


###############################################################################
# 12. Final summary
###############################################################################

summary_lines <- c(
  paste0("Input DCC files: ", length(DCCFiles)),
  paste0("Input PKC files: ", length(PKCFiles)),
  paste0("ROIs after segment QC: ", ncol(target_geomx)),
  paste0("Targets after filtering (including negative controls): ", nrow(target_geomx)),
  paste0("Biological genes exported: ", length(biological_target_names)),
  paste0("Gene detection threshold: ", GENE_DETECTION_THRESHOLD),
  "Normalization: Q3 (75th percentile)"
)

writeLines(
  summary_lines,
  con = file.path(
    OUTPUT_DIR,
    "GeoMx_preprocessing_summary.txt"
  )
)

message("")
message("GeoMx preprocessing completed.")
message("Results written to: ", normalizePath(OUTPUT_DIR, mustWork = FALSE))
dim(target_geomx)
cat(
  readLines(
    file.path(
      OUTPUT_DIR,
      "GeoMx_preprocessing_summary.txt"
    )
  ),
  sep = "\n"
)
