#!/usr/bin/env Rscript

###############################################################################
# 02_CosMx_celltype_annotation.R
#
# CosMx Spatial Molecular Imaging (SMI)
# Obesity-associated colorectal cancer study
#
# Purpose
# -------
# Annotate CosMx cells by label transfer from the colorectal cancer
# single-cell RNA-seq reference (GSE132465).
#
# Workflow
# -------
# 1. Load preprocessed CosMx Seurat object
# 2. Load scRNA-seq reference object
# 3. Transfer major cell-type labels
# 4. Transfer cell-subtype labels
# 5. Add predicted labels and prediction scores to CosMx metadata
# 6. Generate basic UMAP annotation QC plots
# 7. Save annotated CosMx Seurat object
#
# Original analysis strategy
# --------------------------
# FindTransferAnchors(..., normalization.method = "SCT")
# TransferData(..., refdata = reference$Cell_type,    dims = 1:30)
# TransferData(..., refdata = reference$Cell_subtype, dims = 1:30)
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- Sys.getenv(
  "COSMX_PREPROCESSED_RDS",
  unset = "results/CosMx/01_preprocessing/CosMx_preprocessed.rds"
)

REFERENCE_RDS <- Sys.getenv(
  "CRC_REFERENCE_RDS",
  unset = "path/to/gse132465_annotations.rds"
)

OUTPUT_DIR <- Sys.getenv(
  "COSMX_ANNOTATION_OUTPUT",
  unset = "results/CosMx/02_celltype_annotation"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(OUTPUT_DIR, "plots"),
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. Analysis parameters
###############################################################################

# The original label-transfer workflow used dimensions 1:30.
TRANSFER_DIMS <- 1:30

REFERENCE_CELLTYPE_COLUMN <- "Cell_type"
REFERENCE_SUBTYPE_COLUMN  <- "Cell_subtype"


###############################################################################
# 2. Packages
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})


###############################################################################
# 3. Load CosMx query object
###############################################################################

if (!file.exists(INPUT_RDS)) {
  stop(
    "CosMx input RDS was not found:\n",
    INPUT_RDS
  )
}

cosmx <- readRDS(INPUT_RDS)

if (!inherits(cosmx, "Seurat")) {
  stop("INPUT_RDS does not contain a Seurat object.")
}

if (!"SCT" %in% Assays(cosmx)) {
  stop(
    "The CosMx object does not contain an SCT assay. ",
    "Run 01_CosMx_preprocessing.R first."
  )
}

DefaultAssay(cosmx) <- "SCT"

cat(
  "\n============================================================\n",
  "Loaded CosMx query object\n",
  "============================================================\n",
  "Cells: ", ncol(cosmx), "\n",
  "Features: ", nrow(cosmx), "\n",
  sep = ""
)


###############################################################################
# 4. Load scRNA-seq reference
###############################################################################

if (!file.exists(REFERENCE_RDS)) {
  stop(
    "Reference RDS was not found:\n",
    REFERENCE_RDS
  )
}

crc_reference <- readRDS(REFERENCE_RDS)

if (!inherits(crc_reference, "Seurat")) {
  stop("REFERENCE_RDS does not contain a Seurat object.")
}

required_reference_columns <- c(
  REFERENCE_CELLTYPE_COLUMN,
  REFERENCE_SUBTYPE_COLUMN
)

missing_reference_columns <- setdiff(
  required_reference_columns,
  colnames(crc_reference@meta.data)
)

if (length(missing_reference_columns) > 0) {
  stop(
    "Missing reference metadata column(s): ",
    paste(missing_reference_columns, collapse = ", ")
  )
}

if (!"SCT" %in% Assays(crc_reference)) {
  stop(
    "The reference object does not contain an SCT assay, but the original ",
    "label-transfer workflow used normalization.method = 'SCT'."
  )
}

DefaultAssay(crc_reference) <- "SCT"

cat(
  "\n============================================================\n",
  "Loaded GSE132465 reference object\n",
  "============================================================\n",
  "Reference cells: ", ncol(crc_reference), "\n",
  "Major cell types: ",
  length(unique(crc_reference@meta.data[[REFERENCE_CELLTYPE_COLUMN]])),
  "\n",
  "Cell subtypes: ",
  length(unique(crc_reference@meta.data[[REFERENCE_SUBTYPE_COLUMN]])),
  "\n",
  sep = ""
)


###############################################################################
# 5. Major cell-type label transfer
###############################################################################

cat("\nFinding anchors for major cell-type transfer...\n")

anchors_celltype <- FindTransferAnchors(
  reference = crc_reference,
  query = cosmx,
  normalization.method = "SCT",
  dims = TRANSFER_DIMS
)

cat("Transferring major cell-type labels...\n")

prediction_celltype <- TransferData(
  anchorset = anchors_celltype,
  refdata = crc_reference@meta.data[[REFERENCE_CELLTYPE_COLUMN]],
  dims = TRANSFER_DIMS
)

celltype_metadata <- data.frame(
  predicted_cell_type = prediction_celltype$predicted.id,
  prediction_score_cell_type = prediction_celltype$prediction.score.max,
  row.names = rownames(prediction_celltype),
  check.names = FALSE
)

cosmx <- AddMetaData(
  cosmx,
  metadata = celltype_metadata
)


###############################################################################
# 6. Cell-subtype label transfer
###############################################################################

cat("\nFinding anchors for cell-subtype transfer...\n")

anchors_subtype <- FindTransferAnchors(
  reference = crc_reference,
  query = cosmx,
  normalization.method = "SCT",
  dims = TRANSFER_DIMS
)

cat("Transferring cell-subtype labels...\n")

prediction_subtype <- TransferData(
  anchorset = anchors_subtype,
  refdata = crc_reference@meta.data[[REFERENCE_SUBTYPE_COLUMN]],
  dims = TRANSFER_DIMS
)

subtype_metadata <- data.frame(
  predicted_cell_subtype = prediction_subtype$predicted.id,
  prediction_score_cell_subtype = prediction_subtype$prediction.score.max,
  row.names = rownames(prediction_subtype),
  check.names = FALSE
)

cosmx <- AddMetaData(
  cosmx,
  metadata = subtype_metadata
)


###############################################################################
# 7. Annotation summary
###############################################################################

cat("\nMajor cell-type predictions:\n")
print(
  sort(
    table(cosmx$predicted_cell_type),
    decreasing = TRUE
  )
)

cat("\nCell-subtype predictions:\n")
print(
  sort(
    table(cosmx$predicted_cell_subtype),
    decreasing = TRUE
  )
)


celltype_summary <- as.data.frame(
  table(cosmx$predicted_cell_type),
  stringsAsFactors = FALSE
)

colnames(celltype_summary) <- c(
  "predicted_cell_type",
  "cell_count"
)

celltype_summary <- celltype_summary[
  order(
    celltype_summary$cell_count,
    decreasing = TRUE
  ),
]

write.csv(
  celltype_summary,
  file.path(
    OUTPUT_DIR,
    "CosMx_predicted_celltype_counts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


subtype_summary <- as.data.frame(
  table(cosmx$predicted_cell_subtype),
  stringsAsFactors = FALSE
)

colnames(subtype_summary) <- c(
  "predicted_cell_subtype",
  "cell_count"
)

subtype_summary <- subtype_summary[
  order(
    subtype_summary$cell_count,
    decreasing = TRUE
  ),
]

write.csv(
  subtype_summary,
  file.path(
    OUTPUT_DIR,
    "CosMx_predicted_cellsubtype_counts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 8. Basic annotation QC plots
###############################################################################

if ("umap" %in% Reductions(cosmx)) {

  p_celltype <- DimPlot(
    cosmx,
    reduction = "umap",
    group.by = "predicted_cell_type",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("CosMx label transfer: major cell type")

  ggsave(
    file.path(
      OUTPUT_DIR,
      "plots",
      "CosMx_UMAP_predicted_cell_type.png"
    ),
    plot = p_celltype,
    width = 10,
    height = 8,
    dpi = 300
  )


  p_subtype <- DimPlot(
    cosmx,
    reduction = "umap",
    group.by = "predicted_cell_subtype",
    label = FALSE,
    raster = TRUE
  ) +
    ggtitle("CosMx label transfer: cell subtype")

  ggsave(
    file.path(
      OUTPUT_DIR,
      "plots",
      "CosMx_UMAP_predicted_cell_subtype.png"
    ),
    plot = p_subtype,
    width = 12,
    height = 9,
    dpi = 300
  )
}


###############################################################################
# 9. Prediction-score QC
###############################################################################

score_summary <- data.frame(
  annotation = c(
    "Cell_type",
    "Cell_subtype"
  ),
  mean_prediction_score = c(
    mean(
      cosmx$prediction_score_cell_type,
      na.rm = TRUE
    ),
    mean(
      cosmx$prediction_score_cell_subtype,
      na.rm = TRUE
    )
  ),
  median_prediction_score = c(
    median(
      cosmx$prediction_score_cell_type,
      na.rm = TRUE
    ),
    median(
      cosmx$prediction_score_cell_subtype,
      na.rm = TRUE
    )
  )
)

write.csv(
  score_summary,
  file.path(
    OUTPUT_DIR,
    "CosMx_label_transfer_score_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 10. Set major transferred label as identity
###############################################################################

Idents(cosmx) <- "predicted_cell_type"


###############################################################################
# 11. Save annotated CosMx object
###############################################################################

saveRDS(
  cosmx,
  file = file.path(
    OUTPUT_DIR,
    "CosMx_label_transferred.rds"
  )
)

write.csv(
  cosmx@meta.data,
  file.path(
    OUTPUT_DIR,
    "CosMx_label_transferred_metadata.csv"
  ),
  row.names = TRUE,
  quote = FALSE
)


###############################################################################
# 12. Session information
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)


###############################################################################
# 13. Final summary
###############################################################################

cat(
  "\n============================================================\n",
  "CosMx label transfer completed\n",
  "============================================================\n",
  "Reference: GSE132465\n",
  "Normalization method: SCT\n",
  "Transfer dimensions: 1:30\n",
  "Major label column: predicted_cell_type\n",
  "Subtype label column: predicted_cell_subtype\n",
  "Output object: ",
  file.path(
    OUTPUT_DIR,
    "CosMx_label_transferred.rds"
  ),
  "\n",
  "============================================================\n",
  sep = ""
)
