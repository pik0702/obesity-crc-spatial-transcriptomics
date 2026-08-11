#!/usr/bin/env Rscript

###############################################################################
# 01_CosMx_preprocessing.R
#
# CosMx Spatial Molecular Imaging (SMI)
# Obesity-associated colorectal cancer study
#
#
# Workflow
# -------
# 1. Load merged Seurat object
# 2. Perform basic cell-level QC
# 3. Remove CosMx control targets
# 4. SCTransform normalization
# 5. PCA
# 6. Harmony integration by patient
# 7. UMAP
# 8. Neighbor graph construction and clustering
# 9. Generate basic UMAP QC plots
# 10. Save the processed Seurat object
#
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- Sys.getenv(
  "COSMX_MERGED_RDS",
  unset = "path/to/CosMx_merged_object.rds"
)

OUTPUT_DIR <- Sys.getenv(
  "COSMX_PREPROCESSING_OUTPUT",
  unset = "results/CosMx/01_preprocessing"
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

MIN_NCOUNT_RNA <- 20

SCT_CLIP_RANGE <- c(-10, 10)

N_PCS <- 27

HARMONY_VARIABLE <- "patient"

CLUSTER_RESOLUTION <- 0.3


###############################################################################
# 2. Packages
###############################################################################

required_packages <- c(
  "Seurat",
  "harmony",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})


###############################################################################
# 3. Helper functions
###############################################################################

remove_cosmx_control_features <- function(object, assay = "RNA") {

  if (!assay %in% Assays(object)) {
    stop("Assay '", assay, "' was not found in the Seurat object.")
  }

  feature_names <- rownames(object[[assay]])

  control_features <- feature_names[
    grepl(
      "^(SystemControl|Negative)",
      feature_names,
      ignore.case = FALSE
    )
  ]

  if (length(control_features) == 0) {
    message("No SystemControl/Negative features detected.")
    return(object)
  }

  biological_features <- setdiff(
    feature_names,
    control_features
  )

  message(
    "Removing ",
    length(control_features),
    " CosMx control features."
  )

  subset(
    object,
    features = biological_features
  )
}


save_plot <- function(plot_object, filename, width = 8, height = 6) {

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "plots",
      filename
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}


###############################################################################
# 4. Load merged CosMx Seurat object
###############################################################################

if (!file.exists(INPUT_RDS)) {
  stop(
    "Input RDS was not found:\n",
    INPUT_RDS
  )
}

cosmx <- readRDS(INPUT_RDS)

if (!inherits(cosmx, "Seurat")) {
  stop("The input RDS does not contain a Seurat object.")
}

cat(
  "\n============================================================\n",
  "Loaded merged CosMx Seurat object\n",
  "============================================================\n",
  "Cells: ", ncol(cosmx), "\n",
  "Assays: ", paste(Assays(cosmx), collapse = ", "), "\n",
  sep = ""
)


###############################################################################
# 5. Basic QC
###############################################################################

if (!"RNA" %in% Assays(cosmx)) {
  stop("RNA assay is required for preprocessing.")
}

DefaultAssay(cosmx) <- "RNA"

if (!"nCount_RNA" %in% colnames(cosmx@meta.data)) {
  stop("Metadata column 'nCount_RNA' was not found.")
}

if (!HARMONY_VARIABLE %in% colnames(cosmx@meta.data)) {
  stop(
    "Metadata column '",
    HARMONY_VARIABLE,
    "' required for Harmony was not found."
  )
}

qc_before <- data.frame(
  cells_before = ncol(cosmx),
  median_nCount_RNA = median(
    cosmx$nCount_RNA,
    na.rm = TRUE
  ),
  median_nFeature_RNA = if (
    "nFeature_RNA" %in% colnames(cosmx@meta.data)
  ) {
    median(
      cosmx$nFeature_RNA,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
)

# Historical cell-level QC threshold used in the CosMx workflow
cosmx <- subset(
  cosmx,
  subset = nCount_RNA > MIN_NCOUNT_RNA
)

qc_after <- data.frame(
  cells_after = ncol(cosmx),
  cells_removed = qc_before$cells_before - ncol(cosmx)
)

qc_summary <- cbind(
  qc_before,
  qc_after
)

write.csv(
  qc_summary,
  file.path(
    OUTPUT_DIR,
    "CosMx_basic_QC_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("\nBasic QC summary:\n")
print(qc_summary)


###############################################################################
# 6. Remove CosMx control targets
###############################################################################

features_before <- nrow(cosmx[["RNA"]])

cosmx <- remove_cosmx_control_features(
  cosmx,
  assay = "RNA"
)

features_after <- nrow(cosmx[["RNA"]])

control_summary <- data.frame(
  features_before = features_before,
  features_after = features_after,
  features_removed = features_before - features_after
)

write.csv(
  control_summary,
  file.path(
    OUTPUT_DIR,
    "CosMx_control_feature_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 7. SCTransform
###############################################################################

DefaultAssay(cosmx) <- "RNA"

cosmx <- SCTransform(
  cosmx,
  assay = "RNA",
  clip.range = SCT_CLIP_RANGE,
  verbose = FALSE
)


###############################################################################
# 8. PCA
###############################################################################

cosmx <- RunPCA(
  cosmx,
  assay = "SCT",
  npcs = N_PCS,
  verbose = FALSE
)


###############################################################################
# 9. Harmony integration
#
# Final workflow:
#   group.by.vars = "patient"
###############################################################################

cosmx <- harmony::RunHarmony(
  object = cosmx,
  group.by.vars = HARMONY_VARIABLE,
  reduction = "pca"
)


###############################################################################
# 10. UMAP
#
# Final workflow:
#   reduction = "harmony"
#   dims = 1:27
###############################################################################

cosmx <- RunUMAP(
  cosmx,
  reduction = "harmony",
  dims = seq_len(N_PCS),
  verbose = FALSE
)


###############################################################################
# 11. Neighbor graph and clustering
#
# Final workflow:
#   reduction = "harmony"
#   dims = 1:27
#   resolution = 0.3
###############################################################################

cosmx <- FindNeighbors(
  cosmx,
  reduction = "harmony",
  dims = seq_len(N_PCS),
  verbose = FALSE
)

cosmx <- FindClusters(
  cosmx,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)


###############################################################################
# 12. Basic UMAP QC plots
###############################################################################

p_cluster <- DimPlot(
  cosmx,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("CosMx UMAP by cluster")

save_plot(
  p_cluster,
  "CosMx_UMAP_by_cluster.png"
)


p_patient <- DimPlot(
  cosmx,
  reduction = "umap",
  group.by = "patient",
  label = FALSE
) +
  ggtitle("CosMx UMAP by patient")

save_plot(
  p_patient,
  "CosMx_UMAP_by_patient.png",
  width = 10,
  height = 8
)


if ("obese_status" %in% colnames(cosmx@meta.data)) {

  p_obesity <- DimPlot(
    cosmx,
    reduction = "umap",
    group.by = "obese_status",
    label = FALSE
  ) +
    ggtitle("CosMx UMAP by obesity status")

  save_plot(
    p_obesity,
    "CosMx_UMAP_by_obesity_status.png"
  )
}


###############################################################################
# 13. Save processed object and metadata
###############################################################################

saveRDS(
  cosmx,
  file = file.path(
    OUTPUT_DIR,
    "CosMx_preprocessed.rds"
  )
)

write.csv(
  cosmx@meta.data,
  file.path(
    OUTPUT_DIR,
    "CosMx_preprocessed_cell_metadata.csv"
  ),
  row.names = TRUE,
  quote = FALSE
)


###############################################################################
# 14. Session information
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)


###############################################################################
# 15. Final summary
###############################################################################

cat(
  "\n============================================================\n",
  "CosMx preprocessing completed\n",
  "============================================================\n",
  "Cells: ", ncol(cosmx), "\n",
  "RNA features: ", nrow(cosmx[["RNA"]]), "\n",
  "SCT clip range: -10 to 10\n",
  "PCA dimensions: 1-", N_PCS, "\n",
  "Harmony variable: ", HARMONY_VARIABLE, "\n",
  "UMAP/neighbor reduction: harmony\n",
  "Clustering resolution: ", CLUSTER_RESOLUTION, "\n",
  "Output: ", normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  ), "\n",
  "============================================================\n",
  sep = ""
)
