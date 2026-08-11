#!/usr/bin/env Rscript

###############################################################################
# GeoMx DSP visualization
#
# Study:
#   Obesity-associated colorectal cancer spatial transcriptomics
#
#
# Input:
#   results/GeoMx/01_preprocessing_QC_normalization/
#     GeoMx_preprocessed_Q3_normalized.rds
#
# Output:
#   results/GeoMx/06_visualization/
#
# Notes:
#   - Figure 1C follows the original marker-bubble workflow:
#       log2(Q3 + 1) expression
#       mean expression by ROI group
#       gene-wise normalized bubble size
#
#   - Figure 1D follows the original UMAP workflow using log2(Q3)
#     expression and a fixed random seed.
#
#   - Figure 1E follows the original heatmap workflow:
#       log2(Q3) expression
#       coefficient of variation (CV) per gene
#       top 20% most variable genes
#       row scaling and correlation-based hierarchical clustering
#
###############################################################################


###############################################################################
# 0. Configuration
###############################################################################

INPUT_RDS <- Sys.getenv(
  "GEOMX_PREPROCESSED_RDS",
  unset = file.path(
    "results",
    "GeoMx",
    "01_preprocessing_QC_normalization",
    "GeoMx_preprocessed_Q3_normalized.rds"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "GEOMX_VIS_OUTPUT_DIR",
  unset = file.path(
    "results",
    "GeoMx",
    "06_visualization"
  )
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. Required packages
###############################################################################

required_packages <- c(
  "Biobase",
  "GeomxTools",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "umap",
  "pheatmap"
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
  library(Biobase)
  library(GeomxTools)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(umap)
  library(pheatmap)
})


###############################################################################
# 2. Load preprocessed GeoMx object
###############################################################################

if (!file.exists(INPUT_RDS)) {
  stop(
    "Preprocessed GeoMx object was not found:\n",
    INPUT_RDS,
    "\nRun 01_preprocessing_QC_normalization.R first, or define ",
    "GEOMX_PREPROCESSED_RDS."
  )
}

target_geomx <- readRDS(INPUT_RDS)

if (!"q_norm" %in% Biobase::assayDataElementNames(target_geomx)) {
  stop(
    "The input GeoMx object does not contain the 'q_norm' assay element."
  )
}

q3_matrix <- Biobase::assayDataElement(
  target_geomx,
  elt = "q_norm"
)

roi_metadata <- as.data.frame(
  Biobase::pData(target_geomx),
  check.names = FALSE
)

required_metadata <- c(
  "ROIType",
  "class",
  "Category"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(roi_metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    "Required ROI metadata column(s) are missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

message(
  "Loaded: ",
  nrow(q3_matrix), " targets x ",
  ncol(q3_matrix), " ROIs"
)


###############################################################################
# 3. Standardize display labels
###############################################################################

roi_metadata$ROIType_plot <- dplyr::recode(
  as.character(roi_metadata$ROIType),
  "CD3" = "Immune",
  "Epithelial" = "Epithelial",
  "Adipose" = "Adipose",
  .default = as.character(roi_metadata$ROIType)
)

roi_metadata$Class_plot <- dplyr::recode(
  as.character(roi_metadata$class),
  "normal" = "Normal",
  "Normal" = "Normal",
  "tumor" = "Tumor",
  "Tumor" = "Tumor",
  .default = as.character(roi_metadata$class)
)

# In the original metadata, Category was used to distinguish obesity status.
# "Normal" is converted here to "Non-obese" only for figure labeling.
roi_metadata$Obesity_plot <- dplyr::recode(
  as.character(roi_metadata$Category),
  "Normal" = "Non-obese",
  "Non-obese" = "Non-obese",
  "Nonobese" = "Non-obese",
  "Obese" = "Obese",
  .default = as.character(roi_metadata$Category)
)

rownames(roi_metadata) <- colnames(q3_matrix)


###############################################################################
# 4. Figure 1C - Canonical marker bubble plot
###############################################################################

# Marker order follows the final Figure 1C.
marker_genes <- c(
  # Tumor-enriched epithelial markers
  "MMP7",
  "CLDN1",
  "KRT80",
  
  # Normal epithelial markers
  "MUC2",
  "MUC1",
  "GPA33",
  
  # Shared epithelial markers
  "EPCAM",
  "KRT20",
  
  # Adipose markers
  "FABP4",
  "ADIPOQ",
  
  # Immune markers
  "CD3D",
  "CD3E",
  "IL2RB",
  "CD4",
  "CD8A",
  "CD79A",
  "CD68",
  "CD163",
  "S100A9"
)

marker_group <- c(
  rep("Tumor epithelial", 3),
  rep("Normal epithelial", 3),
  rep("Common epithelial", 2),
  rep("Adipose", 2),
  rep("Immune", 9)
)

missing_markers <- setdiff(
  marker_genes,
  rownames(q3_matrix)
)

if (length(missing_markers) > 0) {
  warning(
    "Marker genes absent from the Q3 matrix and omitted: ",
    paste(missing_markers, collapse = ", ")
  )
}

marker_genes_present <- intersect(
  marker_genes,
  rownames(q3_matrix)
)

if (length(marker_genes_present) == 0) {
  stop("None of the Figure 1C marker genes were found in the Q3 matrix.")
}

marker_group_df <- data.frame(
  Gene = marker_genes,
  MarkerGroup = marker_group,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(Gene %in% marker_genes_present)

# log2(Q3 + 1), matching the final marker-bubble source script.
marker_expr <- log2(
  q3_matrix[
    marker_genes_present,
    ,
    drop = FALSE
  ] + 1
)

marker_long <- as.data.frame(
  marker_expr,
  check.names = FALSE
) |>
  tibble::rownames_to_column("Gene") |>
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "ROI",
    values_to = "Expression"
  ) |>
  dplyr::left_join(
    roi_metadata |>
      tibble::rownames_to_column("ROI") |>
      dplyr::select(
        ROI,
        ROIType_plot,
        Class_plot
      ),
    by = "ROI"
  ) |>
  dplyr::mutate(
    FigureGroup = dplyr::case_when(
      ROIType_plot == "Epithelial" &
        Class_plot == "Tumor" ~ "Epithelial-Tumor",
      
      ROIType_plot == "Epithelial" &
        Class_plot == "Normal" ~ "Epithelial-Normal",
      
      ROIType_plot == "Immune" ~ "Immune",
      ROIType_plot == "Adipose" ~ "Adipose",
      
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(
    !is.na(FigureGroup)
  )

marker_summary <- marker_long |>
  dplyr::group_by(
    Gene,
    FigureGroup
  ) |>
  dplyr::summarise(
    MeanExpression = mean(
      Expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    marker_group_df,
    by = "Gene"
  ) |>
  dplyr::group_by(Gene) |>
  dplyr::mutate(
    min_expr = min(
      MeanExpression,
      na.rm = TRUE
    ),
    max_expr = max(
      MeanExpression,
      na.rm = TRUE
    ),
    NormalizedSize = if (
      max_expr[1] == min_expr[1]
    ) {
      rep(
        5.5,
        dplyr::n()
      )
    } else {
      (
        MeanExpression - min_expr
      ) /
        (
          max_expr - min_expr
        ) * 9 + 1
    }
  ) |>
  dplyr::select(
    -min_expr,
    -max_expr
  ) |>
  dplyr::ungroup()

# ggplot places the first discrete y level at the bottom.
# This order therefore gives the final top-to-bottom display:
# Adipose, Immune, Epithelial(T), Epithelial(N).
marker_summary$FigureGroup <- factor(
  marker_summary$FigureGroup,
  levels = c(
    "Epithelial-Normal",
    "Epithelial-Tumor",
    "Immune",
    "Adipose"
  )
)


marker_summary$Gene <- factor(
  as.character(marker_summary$Gene),
  levels = marker_genes
)
###############################################################################
# Force Figure 1C gene order
###############################################################################

marker_genes <- c(
  "MMP7",
  "CLDN1",
  "KRT80",
  
  "MUC2",
  "MUC1",
  "GPA33",
  
  "EPCAM",
  "KRT20",
  
  "FABP4",
  "ADIPOQ",
  
  "CD3D",
  "CD3E",
  "IL2RB",
  "CD4",
  "CD8A",
  "CD79A",
  "CD68",
  "CD163",
  "S100A9"
)
# Vertical dashed lines separating marker classes.
group_boundaries <- c(
  3.5,   # after KRT80
  6.5,   # after GPA33
  8.5,   # after KRT20
  10.5   # after ADIPOQ
)

fig1c <- ggplot(
  marker_summary,
  aes(
    x = Gene,
    y = FigureGroup,
    size = NormalizedSize,
    color = MeanExpression
  )
) +
  geom_point(alpha = 0.80) +
  
  geom_vline(
    xintercept = c(3.5, 6.5, 8.5, 10.5),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  
  scale_x_discrete(
    limits = marker_genes,
    drop = FALSE
  ) +
  
  scale_size_continuous(
    range = c(1.2, 6.2),
    breaks = c(2.5, 5.0, 7.5, 10.0),
    limits = c(1, 10),
    name = "Normalized\nSize"
  ) +
  
  scale_color_gradient(
    low = "lightblue",
    high = "red",
    name = "Mean\nExpression"
  ) +
  
  scale_y_discrete(
    labels = c(
      "Epithelial-Normal" = "Epithelial(N)",
      "Epithelial-Tumor" = "Epithelial(T)",
      "Immune" = "Immune",
      "Adipose" = "Adipose"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

print(fig1c)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure1C_GeoMx_marker_bubble.pdf"
  ),
  plot = fig1c,
  width = 9.0,
  height = 4.1,
  units = "in"
)

write.csv(
  marker_summary,
  file.path(
    OUTPUT_DIR,
    "Figure1C_marker_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 5. Figure 1D - UMAP of GeoMx ROIs
###############################################################################

if (any(q3_matrix <= 0, na.rm = TRUE)) {
  warning(
    "Non-positive Q3 values were detected. ",
    "For UMAP, log2(Q3 + 1) will be used instead of log2(Q3)."
  )
  
  umap_input <- log2(
    q3_matrix + 1
  )
  
} else {
  
  # Original workflow used log2(Q3).
  umap_input <- log2(
    q3_matrix
  )
}

set.seed(42)

custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_out <- umap::umap(
  t(umap_input),
  config = custom_umap
)

umap_df <- data.frame(
  ROI = colnames(q3_matrix),
  UMAP1 = umap_out$layout[, 1],
  UMAP2 = umap_out$layout[, 2],
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(
    roi_metadata |>
      tibble::rownames_to_column("ROI") |>
      dplyr::select(
        ROI,
        ROIType_plot,
        Class_plot,
        Obesity_plot
      ),
    by = "ROI"
  )

umap_df$ROIType_plot <- factor(
  umap_df$ROIType_plot,
  levels = c(
    "Epithelial",
    "Immune",
    "Adipose"
  )
)

umap_df$Class_plot <- factor(
  umap_df$Class_plot,
  levels = c(
    "Normal",
    "Tumor"
  )
)

fig1d <- ggplot(
  umap_df,
  aes(
    # Swapped only for display to match the final panel orientation.
    x = UMAP2,
    y = UMAP1,
    color = ROIType_plot,
    shape = Class_plot
  )
) +
  geom_point(
    size = 2.7,
    alpha = 0.85
  ) +
  stat_ellipse(
    aes(
      group = interaction(
        ROIType_plot,
        Class_plot,
        drop = TRUE
      )
    ),
    type = "norm",
    level = 0.80,
    linetype = "dotted",
    linewidth = 0.55,
    color = "black",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Epithelial" = "#619CFF",
      "Immune" = "#00BA38",
      "Adipose" = "#F8766D"
    ),
    name = "ROI type"
  ) +
  scale_shape_manual(
    values = c(
      "Normal" = 16,
      "Tumor" = 17
    ),
    name = "Class"
  ) +
  labs(
    x = "UMAP2",
    y = "UMAP1"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

print(fig1d)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure1D_GeoMx_UMAP.pdf"
  ),
  plot = fig1d,
  width = 6.2,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure1D_GeoMx_UMAP.png"
  ),
  plot = fig1d,
  width = 6.2,
  height = 5.2,
  units = "in",
  dpi = 600
)

write.csv(
  umap_df,
  file.path(
    OUTPUT_DIR,
    "Figure1D_UMAP_coordinates.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 6. Figure 1E - Unsupervised expression heatmap
###############################################################################

if (any(q3_matrix <= 0, na.rm = TRUE)) {
  
  # Avoid undefined log2 values if any non-positive Q3 values remain.
  heatmap_log_expr <- log2(
    q3_matrix + 1
  )
  
} else {
  
  # Original Figure 1E source used log2(Q3).
  heatmap_log_expr <- log2(
    q3_matrix
  )
}

calc_CV <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (
    length(x) < 2 ||
    mean(x) == 0
  ) {
    return(NA_real_)
  }
  
  stats::sd(x) /
    mean(x)
}

CV_dat <- apply(
  heatmap_log_expr,
  1,
  calc_CV
)

CV_cutoff <- stats::quantile(
  CV_dat,
  probs = 0.80,
  na.rm = TRUE
)

GOI <- names(CV_dat)[
  !is.na(CV_dat) &
    CV_dat > CV_cutoff
]

if (length(GOI) == 0) {
  stop(
    "No genes passed the top-20% CV criterion for Figure 1E."
  )
}

heatmap_matrix <- heatmap_log_expr[
  GOI,
  ,
  drop = FALSE
]

# Column annotations in the final figure.
heatmap_annotation <- data.frame(
  `ROI type` = roi_metadata$ROIType_plot,
  Class = roi_metadata$Class_plot,
  Obesity = roi_metadata$Obesity_plot,
  row.names = rownames(roi_metadata),
  check.names = FALSE
)

heatmap_annotation <- heatmap_annotation[
  colnames(heatmap_matrix),
  ,
  drop = FALSE
]

annotation_colors <- list(
  `ROI type` = c(
    "Epithelial" = "#8ECAE6",
    "Immune" = "#F28482",
    "Adipose" = "#F4C95D"
  ),
  Class = c(
    "Normal" = "#BDBDBD",
    "Tumor" = "#E31A1C"
  ),
  Obesity = c(
    "Non-obese" = "#33A6B8",
    "Obese" = "#D8E219"
  )
)

# Remove annotation color levels that are absent in the current dataset.
for (nm in names(annotation_colors)) {
  
  observed <- unique(
    as.character(
      heatmap_annotation[[nm]]
    )
  )
  
  annotation_colors[[nm]] <- annotation_colors[[nm]][
    names(annotation_colors[[nm]]) %in% observed
  ]
}

# Save selected variable genes and their CV values.
variable_gene_table <- data.frame(
  Gene = GOI,
  CV = CV_dat[GOI],
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(
    dplyr::desc(CV)
  )

write.csv(
  variable_gene_table,
  file.path(
    OUTPUT_DIR,
    "Figure1E_top20pct_CV_genes.csv"
  ),
  row.names = FALSE
)

# PDF
pdf(
  file.path(
    OUTPUT_DIR,
    "Figure1E_GeoMx_expression_heatmap.pdf"
  ),
  width = 7.0,
  height = 8.0
)

pheatmap::pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  breaks = seq(
    -3,
    3,
    0.05
  ),
  color = grDevices::colorRampPalette(
    c(
      "purple3",
      "black",
      "yellow2"
    )
  )(120),
  annotation_col = heatmap_annotation,
  annotation_colors = annotation_colors,
  annotation_names_col = FALSE,
  silent = FALSE
)

dev.off()

# PNG
png(
  file.path(
    OUTPUT_DIR,
    "Figure1E_GeoMx_expression_heatmap.png"
  ),
  width = 4200,
  height = 4800,
  res = 600
)

pheatmap::pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  breaks = seq(
    -3,
    3,
    0.05
  ),
  color = grDevices::colorRampPalette(
    c(
      "purple3",
      "black",
      "yellow2"
    )
  )(120),
  annotation_col = heatmap_annotation,
  annotation_colors = annotation_colors,
  annotation_names_col = FALSE,
  silent = FALSE
)

dev.off()


###############################################################################
# 7. Save figure-generation metadata
###############################################################################

summary_lines <- c(
  paste0(
    "Input object: ",
    INPUT_RDS
  ),
  paste0(
    "Targets: ",
    nrow(q3_matrix)
  ),
  paste0(
    "ROIs: ",
    ncol(q3_matrix)
  ),
  paste0(
    "Figure 1C markers requested: ",
    length(marker_genes)
  ),
  paste0(
    "Figure 1C markers present: ",
    length(marker_genes_present)
  ),
  paste0(
    "Figure 1E CV cutoff (80th percentile): ",
    signif(CV_cutoff, 6)
  ),
  paste0(
    "Figure 1E genes retained: ",
    length(GOI)
  ),
  "UMAP random seed: 42"
)

writeLines(
  summary_lines,
  con = file.path(
    OUTPUT_DIR,
    "Figure1_CDE_generation_summary.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

message("")
message("GeoMx Figure 1C-1E visualization completed.")
message(
  "Results written to: ",
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  )
)
