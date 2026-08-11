#!/usr/bin/env Rscript

###############################################################################
# 02_differential_expression_early_risk.R
#
# GeoMx DSP obesity-associated CRC study
#
# Input
# -----
# Output from:
#   01_preprocessing_QC_normalization.R
#
# Main workflow
# -------------
# 1. Load Q3-normalized GeoMx object
# 2. Use log2(q_norm) expression
# 3. Define epithelial groups:
#      NO- = normal / non-obese
#      NO+ = normal / obese
#      TO- = tumor / non-obese
#      TO+ = tumor / obese
# 4. Run genome-wide differential-expression statistics
# 5. Reconstruct early-risk candidate logic
# 6. Export the final 30-gene panel
# 7. Recompute the Figure 2B fold-change matrix
# 8. Reproduce the final 30-gene Table S2-2 statistics
# 9. Draw four-group epithelial boxplots (NO-, NO+, TO-, TO+)
#
# Notes
# -----
# - The current public workflow uses all 288 QC-passed ROIs from script 01.
# - The final 30-gene panel is explicitly recorded because biological relevance
#   was part of the final study-level selection.
# - Table S2-2 is calculated with the same historical code structure used for
#   the manuscript figure/statistics.
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- Sys.getenv(
  "GEOMX_PREPROCESSED_RDS",
  unset = paste0(
    "C:/Users/user/Documents/results/GeoMx/",
    "01_preprocessing_QC_normalization/",
    "GeoMx_preprocessed_Q3_normalized.rds"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "GEOMX_FIGURE2_OUTPUT",
  unset = paste0(
    "C:/Users/user/Documents/results/GeoMx/",
    "02_differential_expression_early_risk"
  )
)

RUN_GENOME_WIDE_DGE <- TRUE
MAKE_FIGURE2B_BUBBLE <- TRUE
MAKE_INDIVIDUAL_BOXPLOTS <- TRUE
MAKE_COMBINED_BOXPLOT_PANEL <- TRUE

P_VALUE_CUTOFF <- 0.05
FC_UP <- 1.20
FC_DOWN <- 1 / FC_UP

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    OUTPUT_DIR,
    "boxplots_30genes"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "Biobase",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggpubr",
  "car"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (
  length(
    missing_packages
  ) > 0
) {
  stop(
    "Missing R package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall them first."
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(car)
})


###############################################################################
# 2. Study definitions
###############################################################################

FINAL_30_GENES <- c(
  "NES",
  "PFDN4",
  "CMSS1",
  "DACH1",
  "BGN",
  "CCL24",
  "PTPN13",
  "FGFRL1",
  "KRT6A",
  "MYL9",
  "RAB22A",
  "MCL1",
  "AGR2",
  "FOSL2",
  "C12orf57",
  "NFKBIA",
  "COX5A",
  "SCAMP2",
  "SPINK1",
  "ZFP36",
  "ARL14",
  "BCL10",
  "SYNPO",
  "REG4",
  "AGR3",
  "TMPRSS2",
  "SPON1",
  "MUC2",
  "NANS",
  "SERF2"
)

# Genes displayed as four-group boxplots in the supplied Figure 2 layout.
FIGURE2_BOXPLOT_GENES <- c(
  "NES",
  "CMSS1",
  "PTPN13",
  "FGFRL1",
  "MCL1",
  "SPINK1",
  "DACH1",
  "PFDN4",
  "RAB22A",
  "BGN",
  "CCL24",
  "MYL9",
  "AGR2",
  "MUC2",
  "NANS",
  "COX5A",
  "ARL14",
  "BCL10"
)

GROUP_LEVELS <- c(
  "NO-",
  "NO+",
  "TO-",
  "TO+"
)

COMPARISONS <- list(
  c(
    "NO-",
    "NO+"
  ),
  c(
    "NO-",
    "TO-"
  ),
  c(
    "NO+",
    "TO+"
  ),
  c(
    "TO-",
    "TO+"
  )
)

COMPARISON_NAMES <- c(
  "NO-_vs_NO+",
  "NO-_vs_TO-",
  "NO+_vs_TO+",
  "TO-_vs_TO+"
)


###############################################################################
# 3. Helper functions
###############################################################################

find_column <- function(
    df,
    candidates,
    label
) {

  available <- names(
    df
  )

  index <- match(
    tolower(
      candidates
    ),
    tolower(
      available
    ),
    nomatch = 0
  )

  index <- index[
    index > 0
  ]

  if (
    length(
      index
    ) == 0
  ) {
    stop(
      "Could not identify metadata column for ",
      label,
      ".\nAvailable columns: ",
      paste(
        available,
        collapse = ", "
      )
    )
  }

  available[
    index[1]
  ]
}


standardize_roi_type <- function(
    x
) {

  value <- tolower(
    trimws(
      as.character(
        x
      )
    )
  )

  dplyr::case_when(
    value %in% c(
      "epithelial",
      "epithelium"
    ) ~ "Epithelial",

    value %in% c(
      "cd3",
      "immune"
    ) ~ "CD3",

    value %in% c(
      "adipose",
      "adipose tissue"
    ) ~ "Adipose",

    TRUE ~ as.character(
      x
    )
  )
}


standardize_tissue_class <- function(
    x
) {

  value <- tolower(
    trimws(
      as.character(
        x
      )
    )
  )

  dplyr::case_when(
    value == "normal" ~ "Normal",

    value %in% c(
      "tumor",
      "tumour"
    ) ~ "Tumor",

    TRUE ~ as.character(
      x
    )
  )
}


standardize_obesity <- function(
    x
) {

  value <- tolower(
    trimws(
      as.character(
        x
      )
    )
  )

  dplyr::case_when(
    value %in% c(
      "normal",
      "non-obese",
      "nonobese",
      "non obese"
    ) ~ "Non-obese",

    value %in% c(
      "obese",
      "obesity"
    ) ~ "Obese",

    TRUE ~ as.character(
      x
    )
  )
}


safe_shapiro <- function(
    x
) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (
    length(
      x
    ) < 3 ||
      length(
        unique(
          x
        )
      ) < 2
  ) {
    return(
      0.01
    )
  }

  result <- tryCatch(
    stats::shapiro.test(
      x
    )$p.value,
    error = function(
        e
    ) {
      NA_real_
    }
  )

  if (
    is.na(
      result
    )
  ) {
    result <- 0.01
  }

  result
}


safe_levene <- function(
    x,
    y
) {

  values <- c(
    x,
    y
  )

  groups <- factor(
    c(
      rep(
        "A",
        length(
          x
        )
      ),
      rep(
        "B",
        length(
          y
        )
      )
    )
  )

  keep <- is.finite(
    values
  )

  values <- values[
    keep
  ]

  groups <- droplevels(
    groups[
      keep
    ]
  )

  result <- tryCatch(
    {
      fit <- car::leveneTest(
        values,
        groups,
        center = median
      )

      as.numeric(
        fit[
          1,
          "Pr(>F)"
        ]
      )
    },
    error = function(
        e
    ) {
      NA_real_
    }
  )

  if (
    is.na(
      result
    )
  ) {
    result <- 0.01
  }

  result
}


run_one_gene_test <- function(
    x,
    y
) {

  x <- as.numeric(
    x
  )

  y <- as.numeric(
    y
  )

  x <- x[
    is.finite(
      x
    )
  ]

  y <- y[
    is.finite(
      y
    )
  ]

  shapiro_x <- safe_shapiro(
    x
  )

  shapiro_y <- safe_shapiro(
    y
  )

  levene_p <- safe_levene(
    x,
    y
  )

  both_normal <-
    shapiro_x > 0.05 &&
    shapiro_y > 0.05

  if (
    both_normal &&
      levene_p > 0.05
  ) {

    test_name <- "Student_t"

    p_value <- tryCatch(
      stats::t.test(
        x,
        y,
        var.equal = TRUE
      )$p.value,
      error = function(
          e
      ) {
        NA_real_
      }
    )

  } else if (
    both_normal &&
      levene_p <= 0.05
  ) {

    test_name <- "Welch_t"

    p_value <- tryCatch(
      stats::t.test(
        x,
        y,
        var.equal = FALSE
      )$p.value,
      error = function(
          e
      ) {
        NA_real_
      }
    )

  } else {

    test_name <- "Wilcoxon"

    p_value <- tryCatch(
      suppressWarnings(
        stats::wilcox.test(
          x,
          y,
          exact = FALSE
        )$p.value
      ),
      error = function(
          e
      ) {
        NA_real_
      }
    )
  }

  median_x <- stats::median(
    x,
    na.rm = TRUE
  )

  median_y <- stats::median(
    y,
    na.rm = TRUE
  )

  log2_fc <- median_y -
    median_x

  fold_change <- 2^log2_fc

  c(
    Shapiro_Group1 = shapiro_x,
    Shapiro_Group2 = shapiro_y,
    Levene_p = levene_p,
    Test = test_name,
    Median_log2_Group1 = median_x,
    Median_log2_Group2 = median_y,
    FoldChange = fold_change,
    log2FC = log2_fc,
    p_value = p_value
  )
}


run_genomewide_comparison <- function(
    expression_matrix,
    group_vector,
    group1,
    group2,
    comparison_name
) {

  sample1 <- names(
    group_vector
  )[
    group_vector ==
      group1
  ]

  sample2 <- names(
    group_vector
  )[
    group_vector ==
      group2
  ]

  if (
    length(
      sample1
    ) == 0 ||
      length(
        sample2
      ) == 0
  ) {
    stop(
      "No samples found for comparison ",
      comparison_name
    )
  }

  result_list <- vector(
    "list",
    nrow(
      expression_matrix
    )
  )

  gene_names <- rownames(
    expression_matrix
  )

  for (
    i in seq_len(
      nrow(
        expression_matrix
      )
    )
  ) {

    gene <- gene_names[i]

    test_result <- run_one_gene_test(
      expression_matrix[
        i,
        sample1,
        drop = TRUE
      ],
      expression_matrix[
        i,
        sample2,
        drop = TRUE
      ]
    )

    result_list[i] <- list(
      data.frame(
        Gene = gene,
        Comparison = comparison_name,
        Group1 = group1,
        Group2 = group2,
        Shapiro_Group1 = as.numeric(
          test_result[
            "Shapiro_Group1"
          ]
        ),
        Shapiro_Group2 = as.numeric(
          test_result[
            "Shapiro_Group2"
          ]
        ),
        Levene_p = as.numeric(
          test_result[
            "Levene_p"
          ]
        ),
        Test = as.character(
          test_result[
            "Test"
          ]
        ),
        Median_log2_Group1 = as.numeric(
          test_result[
            "Median_log2_Group1"
          ]
        ),
        Median_log2_Group2 = as.numeric(
          test_result[
            "Median_log2_Group2"
          ]
        ),
        FoldChange = as.numeric(
          test_result[
            "FoldChange"
          ]
        ),
        log2FC = as.numeric(
          test_result[
            "log2FC"
          ]
        ),
        p_value = as.numeric(
          test_result[
            "p_value"
          ]
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  result <- dplyr::bind_rows(
    result_list
  )

  result$FDR <- stats::p.adjust(
    result$p_value,
    method = "BH"
  )

  result
}


direction_class <- function(
    fold_change
) {

  dplyr::case_when(
    fold_change >= FC_UP ~ "Increased",
    fold_change <= FC_DOWN ~ "Decreased",
    TRUE ~ "No_change"
  )
}


make_four_group_boxplot <- function(
    gene,
    epithelial_plot_data
) {

  gene_df <- data.frame(
    Sample = epithelial_plot_data$Sample,
    Group = factor(
      epithelial_plot_data$Group,
      levels = GROUP_LEVELS
    ),
    Expression = as.numeric(
      getElement(
        epithelial_plot_data,
        gene
      )
    ),
    stringsAsFactors = FALSE
  )

  ggplot(
    gene_df,
    aes(
      x = Group,
      y = Expression,
      fill = Group
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65,
      linewidth = 0.55
    ) +
    geom_jitter(
      width = 0.16,
      height = 0,
      size = 1.25,
      alpha = 0.65
    ) +
    ggpubr::stat_compare_means(
      comparisons = COMPARISONS,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE,
      size = 4,
      step.increase = 0.08
    ) +
    scale_fill_manual(
      values = c(
        "NO-" = "#BDBDBD",
        "NO+" = "#FFF4B8",
        "TO-" = "#F3A08E",
        "TO+" = "#F6C4C9"
      ),
      drop = FALSE
    ) +
    labs(
      title = gene,
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.text.x = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.y = element_text(
        size = 9
      ),
      axis.line = element_line(
        linewidth = 0.5
      ),
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    ) +
    coord_cartesian(
      clip = "off"
    )
}


###############################################################################
# 4. Load 01 preprocessing output
###############################################################################

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input RDS not found:\n",
    INPUT_RDS
  )
}

cat(
  "\nLoading GeoMx preprocessing output...\n"
)

geomx <- readRDS(
  INPUT_RDS
)

q_norm_all <- Biobase::assayDataElement(
  geomx,
  elt = "q_norm"
)

phenodata <- as.data.frame(
  Biobase::pData(
    geomx
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

featuredata <- as.data.frame(
  Biobase::fData(
    geomx
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Loaded targets x ROIs: ",
  nrow(
    q_norm_all
  ),
  " x ",
  ncol(
    q_norm_all
  ),
  "\n",
  sep = ""
)


###############################################################################
# 5. Retain biological genes and create log2(Q3)
###############################################################################

biological_keep <- rep(
  TRUE,
  nrow(
    q_norm_all
  )
)

if (
  "CodeClass" %in%
    names(
      featuredata
    )
) {

  biological_keep <- as.character(
    featuredata$CodeClass
  ) != "Negative"
}

gene_names <- if (
  "TargetName" %in%
    names(
      featuredata
    )
) {
  as.character(
    featuredata$TargetName
  )
} else {
  rownames(
    q_norm_all
  )
}

q_norm <- q_norm_all[
  biological_keep,
  ,
  drop = FALSE
]

rownames(
  q_norm
) <- gene_names[
  biological_keep
]

if (
  anyDuplicated(
    rownames(
      q_norm
    )
  ) > 0
) {
  stop(
    "Duplicated biological TargetName values were detected."
  )
}

if (
  any(
    q_norm <= 0,
    na.rm = TRUE
  )
) {
  stop(
    "Non-positive q_norm values detected; log2(q_norm) cannot be calculated."
  )
}

log2_q3 <- log2(
  q_norm
)

cat(
  "Biological genes retained: ",
  nrow(
    log2_q3
  ),
  "\n",
  sep = ""
)


###############################################################################
# 6. Build final epithelial metadata
###############################################################################

if (
  !all(
    colnames(
      log2_q3
    ) %in%
      rownames(
        phenodata
      )
  )
) {
  stop(
    "Expression matrix columns do not match pData row names."
  )
}

phenodata <- phenodata[
  colnames(
    log2_q3
  ),
  ,
  drop = FALSE
]

roi_col <- find_column(
  phenodata,
  c(
    "ROIType",
    "ROI_Type",
    "ROI Type",
    "segment"
  ),
  "ROI type"
)

class_col <- find_column(
  phenodata,
  c(
    "class",
    "Class"
  ),
  "normal/tumor class"
)

obesity_col <- find_column(
  phenodata,
  c(
    "Category",
    "category",
    "Obesity"
  ),
  "obesity category"
)

metadata <- data.frame(
  Sample = rownames(
    phenodata
  ),
  ROIType = standardize_roi_type(
    getElement(
      phenodata,
      roi_col
    )
  ),
  TissueClass = standardize_tissue_class(
    getElement(
      phenodata,
      class_col
    )
  ),
  Obesity = standardize_obesity(
    getElement(
      phenodata,
      obesity_col
    )
  ),
  stringsAsFactors = FALSE
)

metadata$Group <- dplyr::case_when(
  metadata$ROIType == "Epithelial" &
    metadata$TissueClass == "Normal" &
    metadata$Obesity == "Non-obese" ~ "NO-",

  metadata$ROIType == "Epithelial" &
    metadata$TissueClass == "Normal" &
    metadata$Obesity == "Obese" ~ "NO+",

  metadata$ROIType == "Epithelial" &
    metadata$TissueClass == "Tumor" &
    metadata$Obesity == "Non-obese" ~ "TO-",

  metadata$ROIType == "Epithelial" &
    metadata$TissueClass == "Tumor" &
    metadata$Obesity == "Obese" ~ "TO+",

  TRUE ~ NA_character_
)

epithelial_metadata <- metadata[
  metadata$ROIType ==
    "Epithelial" &
    !is.na(
      metadata$Group
    ),
  ,
  drop = FALSE
]

cat(
  "\nEpithelial group counts:\n"
)

print(
  table(
    epithelial_metadata$Group
  )
)

expected_counts <- c(
  "NO-" = 16,
  "NO+" = 24,
  "TO-" = 33,
  "TO+" = 42
)

observed_counts <- table(
  factor(
    epithelial_metadata$Group,
    levels = names(
      expected_counts
    )
  )
)

if (
  !all(
    as.integer(
      observed_counts
    ) ==
      as.integer(
        expected_counts
      )
  )
) {
  warning(
    "Epithelial group counts differ from the final Figure 2 dataset ",
    "(expected 16 / 24 / 33 / 42)."
  )
}

epithelial_samples <- epithelial_metadata$Sample

epithelial_log2_q3 <- log2_q3[
  ,
  epithelial_samples,
  drop = FALSE
]

group_vector <- epithelial_metadata$Group
names(
  group_vector
) <- epithelial_metadata$Sample


###############################################################################
# 7. Genome-wide DGE
#
# Historical statistical family:
#   - Shapiro-Wilk in each group
#   - Levene variance test
#   - Student t-test if both groups are normal and variances are equal
#   - Welch t-test if both groups are normal and variances are unequal
#   - Wilcoxon rank-sum test otherwise
#
# Fold change is calculated on the Q3 scale as:
#   2^(median(log2Q3_Group2) - median(log2Q3_Group1))
###############################################################################

dge_all <- NULL

if (
  RUN_GENOME_WIDE_DGE
) {

  cat(
    "\nRunning genome-wide epithelial DGE...\n"
  )

  dge_list <- vector(
    "list",
    length(
      COMPARISONS
    )
  )

  for (
    i in seq_along(
      COMPARISONS
    )
  ) {

    comparison <- getElement(
      COMPARISONS,
      i
    )

    comparison_name <- COMPARISON_NAMES[i]

    cat(
      "  ",
      comparison_name,
      " ...\n",
      sep = ""
    )

    dge_list[i] <- list(
      run_genomewide_comparison(
        expression_matrix = epithelial_log2_q3,
        group_vector = group_vector,
        group1 = comparison[1],
        group2 = comparison[2],
        comparison_name = comparison_name
      )
    )
  }

  dge_all <- dplyr::bind_rows(
    dge_list
  )

  write.csv(
    dge_all,
    file.path(
      OUTPUT_DIR,
      "Figure2_genomewide_DGE_all_comparisons.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  dge_summary <- dge_all |>
    dplyr::group_by(
      Comparison
    ) |>
    dplyr::summarise(
      N_genes = dplyr::n(),
      Student_t = sum(
        Test ==
          "Student_t",
        na.rm = TRUE
      ),
      Welch_t = sum(
        Test ==
          "Welch_t",
        na.rm = TRUE
      ),
      Wilcoxon = sum(
        Test ==
          "Wilcoxon",
        na.rm = TRUE
      ),
      P005 = sum(
        p_value <
          P_VALUE_CUTOFF,
        na.rm = TRUE
      ),
      FDR005 = sum(
        FDR <
          P_VALUE_CUTOFF,
        na.rm = TRUE
      ),
      P005_FC1.2 = sum(
        p_value <
          P_VALUE_CUTOFF &
          (
            FoldChange >=
              FC_UP |
              FoldChange <=
              FC_DOWN
          ),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    dge_summary,
    file.path(
      OUTPUT_DIR,
      "Figure2_genomewide_DGE_summary.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  cat(
    "\nDGE summary:\n"
  )

  print(
    dge_summary
  )
}


###############################################################################
# 8. Early-risk candidate reconstruction
#
# Case 1:
#   NO- -> NO+ and NO- -> TO- are both significant and direction-consistent.
#
# Case 2:
#   NO- -> NO+ and NO+ -> TO+ are both significant and direction-consistent.
#
# Final panel:
#   30 study-selected genes retained after consistency/biological review.
###############################################################################

early_risk_status <- NULL

if (
  !is.null(
    dge_all
  )
) {

  early <- dge_all |>
    dplyr::filter(
      Comparison ==
        "NO-_vs_NO+"
    ) |>
    dplyr::select(
      Gene,
      early_FC = FoldChange,
      early_p = p_value,
      early_FDR = FDR
    )

  set1 <- dge_all |>
    dplyr::filter(
      Comparison ==
        "NO-_vs_TO-"
    ) |>
    dplyr::select(
      Gene,
      set1_FC = FoldChange,
      set1_p = p_value,
      set1_FDR = FDR
    )

  set2 <- dge_all |>
    dplyr::filter(
      Comparison ==
        "NO+_vs_TO+"
    ) |>
    dplyr::select(
      Gene,
      set2_FC = FoldChange,
      set2_p = p_value,
      set2_FDR = FDR
    )

  late <- dge_all |>
    dplyr::filter(
      Comparison ==
        "TO-_vs_TO+"
    ) |>
    dplyr::select(
      Gene,
      late_FC = FoldChange,
      late_p = p_value,
      late_FDR = FDR
    )

  early_risk_status <- early |>
    dplyr::full_join(
      set1,
      by = "Gene"
    ) |>
    dplyr::full_join(
      set2,
      by = "Gene"
    ) |>
    dplyr::full_join(
      late,
      by = "Gene"
    ) |>
    dplyr::mutate(
      Early_significant =
        early_p <
          P_VALUE_CUTOFF,

      Early_direction =
        direction_class(
          early_FC
        ),

      Set1_direction =
        direction_class(
          set1_FC
        ),

      Set2_direction =
        direction_class(
          set2_FC
        ),

      Case1 =
        Early_significant &
        set1_p <
          P_VALUE_CUTOFF &
        Early_direction %in% c(
          "Increased",
          "Decreased"
        ) &
        Early_direction ==
          Set1_direction,

      Case2 =
        Early_significant &
        set2_p <
          P_VALUE_CUTOFF &
        Early_direction %in% c(
          "Increased",
          "Decreased"
        ) &
        Early_direction ==
          Set2_direction,

      Candidate_union =
        Case1 |
        Case2,

      Final_30 =
        Gene %in%
          FINAL_30_GENES
    )

  write.csv(
    early_risk_status,
    file.path(
      OUTPUT_DIR,
      "Figure2_early_risk_candidate_status.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  final30_status <- early_risk_status |>
    dplyr::filter(
      Final_30
    ) |>
    dplyr::arrange(
      match(
        Gene,
        FINAL_30_GENES
      )
    )

  write.csv(
    final30_status,
    file.path(
      OUTPUT_DIR,
      "Figure2_final_30_gene_status.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  cat(
    "\nEarly-risk reconstruction:\n"
  )

  cat(
    "  Case 1 genes: ",
    sum(
      early_risk_status$Case1,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Case 2 genes: ",
    sum(
      early_risk_status$Case2,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Candidate union: ",
    sum(
      early_risk_status$Candidate_union,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
}


###############################################################################
# 9. Figure 2B-style 30-gene fold-change matrix
#
# This is recomputed consistently from the final 288-ROI object.
# It is intended to reproduce the direction/pattern used in Figure 2B.
###############################################################################

missing_final_genes <- setdiff(
  FINAL_30_GENES,
  rownames(
    epithelial_log2_q3
  )
)

if (
  length(
    missing_final_genes
  ) > 0
) {
  stop(
    "Final panel gene(s) missing from expression matrix: ",
    paste(
      missing_final_genes,
      collapse = ", "
    )
  )
}

figure2b_fc_list <- vector(
  "list",
  length(
    COMPARISONS
  )
)

for (
  i in seq_along(
    COMPARISONS
  )
) {

  comparison <- getElement(
    COMPARISONS,
    i
  )

  group1_samples <- names(
    group_vector
  )[
    group_vector ==
      comparison[1]
  ]

  group2_samples <- names(
    group_vector
  )[
    group_vector ==
      comparison[2]
  ]

  median1 <- apply(
    epithelial_log2_q3[
      FINAL_30_GENES,
      group1_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  )

  median2 <- apply(
    epithelial_log2_q3[
      FINAL_30_GENES,
      group2_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  )

  figure2b_fc_list[i] <- list(
    data.frame(
      Gene = FINAL_30_GENES,
      Comparison = c(
        "NO- -> NO+",
        "NO- -> TO-",
        "NO+ -> TO+",
        "TO- -> TO+"
      )[i],
      log2FC = as.numeric(
        median2 -
          median1
      ),
      FoldChange = as.numeric(
        2^(
          median2 -
            median1
        )
      ),
      stringsAsFactors = FALSE
    )
  )
}

figure2b_long <- dplyr::bind_rows(
  figure2b_fc_list
)

figure2b_long$Direction <- direction_class(
  figure2b_long$FoldChange
)

figure2b_matrix <- figure2b_long |>
  dplyr::select(
    Gene,
    Comparison,
    FoldChange
  ) |>
  tidyr::pivot_wider(
    names_from = Comparison,
    values_from = FoldChange
  ) |>
  dplyr::arrange(
    match(
      Gene,
      FINAL_30_GENES
    )
  )

write.csv(
  figure2b_matrix,
  file.path(
    OUTPUT_DIR,
    "Table_S2-1_recomputed_current_288_ROIs.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  figure2b_long,
  file.path(
    OUTPUT_DIR,
    "Figure2B_fold_change_long.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 10. Figure 2B-style bubble plot
###############################################################################

if (
  MAKE_FIGURE2B_BUBBLE
) {

  bubble_df <- figure2b_long

  bubble_df$Gene <- factor(
    bubble_df$Gene,
    levels = rev(
      FINAL_30_GENES
    )
  )

  bubble_df$Comparison <- factor(
    bubble_df$Comparison,
    levels = c(
      "NO- -> NO+",
      "NO- -> TO-",
      "NO+ -> TO+",
      "TO- -> TO+"
    )
  )

  bubble_plot <- ggplot(
    bubble_df,
    aes(
      x = Comparison,
      y = Gene
    )
  ) +
    geom_point(
      aes(
        fill = log2FC
      ),
      shape = 21,
      size = 4.1,
      color = "grey35",
      stroke = 0.35
    ) +
    scale_fill_gradient2(
      low = "#2C56D8",
      mid = "grey90",
      high = "#E73232",
      midpoint = 0,
      name = "log2 FC"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8
      ),
      axis.ticks = element_blank(),
      legend.position = "right"
    )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "Figure2B_30gene_foldchange_bubble.png"
    ),
    plot = bubble_plot,
    width = 6.5,
    height = 8.5,
    units = "in",
    dpi = 400
  )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "Figure2B_30gene_foldchange_bubble.pdf"
    ),
    plot = bubble_plot,
    width = 6.5,
    height = 8.5,
    units = "in"
  )
}


###############################################################################
# 11. Table S2-2 exact historical calculation
#
# IMPORTANT:
# This block deliberately retains the historical dplyr grouping structure
# used in the final Figure 2 source code. This is the implementation that
# reproduced the manuscript Table S2-2 values from the 01 object.
###############################################################################

epithelial_s22 <- data.frame(
  Sample = epithelial_metadata$Sample,
  Group = as.character(
    epithelial_metadata$Group
  ),
  stringsAsFactors = FALSE
)

for (
  gene in FINAL_30_GENES
) {

  epithelial_s22[
    gene
  ] <- as.numeric(
    epithelial_log2_q3[
      gene,
      epithelial_metadata$Sample,
      drop = TRUE
    ]
  )
}

s22_all_results <- vector(
  "list",
  length(
    COMPARISONS
  )
)

for (
  i in seq_along(
    COMPARISONS
  )
) {

  comparison <- getElement(
    COMPARISONS,
    i
  )

  comparison_name <- COMPARISON_NAMES[i]

  gene_results <- vector(
    "list",
    length(
      FINAL_30_GENES
    )
  )

  for (
    j in seq_along(
      FINAL_30_GENES
    )
  ) {

    gene <- FINAL_30_GENES[j]

    gene_df <- data.frame(
      Group = epithelial_s22$Group,
      Expression = as.numeric(
        getElement(
          epithelial_s22,
          gene
        )
      ),
      stringsAsFactors = FALSE
    )

    gene_df <- gene_df[
      gene_df$Group %in%
        comparison,
      ,
      drop = FALSE
    ]

    test_result <- suppressWarnings(
      stats::wilcox.test(
        Expression ~ Group,
        data = gene_df
      )
    )

    # Preserve the original Figure 2 calculation structure.
    medians <- gene_df |>
      dplyr::group_by(
        Group
      ) |>
      dplyr::summarise(
        med = stats::median(
          Expression,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    fold_change <- medians$med[2] /
      medians$med[1]

    log2_fc <- log2(
      fold_change
    )

    gene_results[j] <- list(
      data.frame(
        Gene = gene,
        Comparison = comparison_name,
        Group1 = comparison[1],
        Group2 = comparison[2],
        Median_Group1 = medians$med[1],
        Median_Group2 = medians$med[2],
        FoldChange = fold_change,
        log2FC = log2_fc,
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
    )
  }

  comparison_df <- dplyr::bind_rows(
    gene_results
  )

  comparison_df$FDR <- stats::p.adjust(
    comparison_df$p_value,
    method = "BH"
  )

  s22_all_results[i] <- list(
    comparison_df
  )
}

table_s22 <- dplyr::bind_rows(
  s22_all_results
)

write.csv(
  table_s22,
  file.path(
    OUTPUT_DIR,
    "Table_S2-2_reproduced.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "\nTable S2-2 rows: ",
  nrow(
    table_s22
  ),
  "\n",
  sep = ""
)


###############################################################################
# 12. Prepare four-group boxplot data
###############################################################################

epithelial_plot_data <- data.frame(
  Sample = epithelial_metadata$Sample,
  Group = as.character(
    epithelial_metadata$Group
  ),
  stringsAsFactors = FALSE
)

for (
  gene in FINAL_30_GENES
) {

  epithelial_plot_data[
    gene
  ] <- as.numeric(
    epithelial_log2_q3[
      gene,
      epithelial_metadata$Sample,
      drop = TRUE
    ]
  )
}

write.csv(
  epithelial_plot_data,
  file.path(
    OUTPUT_DIR,
    "Figure2_boxplot_expression_30genes.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 13. Individual four-group boxplots for all 30 genes
###############################################################################

if (
  MAKE_INDIVIDUAL_BOXPLOTS
) {

  cat(
    "\nDrawing individual four-group boxplots...\n"
  )

  for (
    gene in FINAL_30_GENES
  ) {

    plot_now <- make_four_group_boxplot(
      gene,
      epithelial_plot_data
    )

    ggsave(
      filename = file.path(
        OUTPUT_DIR,
        "boxplots_30genes",
        paste0(
          gene,
          "_NOminus_NOplus_TOminus_TOplus.png"
        )
      ),
      plot = plot_now,
      width = 4,
      height = 4,
      units = "in",
      dpi = 350
    )

    ggsave(
      filename = file.path(
        OUTPUT_DIR,
        "boxplots_30genes",
        paste0(
          gene,
          "_NOminus_NOplus_TOminus_TOplus.pdf"
        )
      ),
      plot = plot_now,
      width = 4,
      height = 4,
      units = "in"
    )
  }
}


###############################################################################
# 14. Combined Figure 2-style boxplot panel
#
# Order follows the supplied Figure 2 example:
#   first 6  : panel C-like genes
#   next 6   : panel D-like genes
#   last 6   : panel E-like genes
###############################################################################

if (
  MAKE_COMBINED_BOXPLOT_PANEL
) {

  cat(
    "Drawing combined Figure 2-style boxplot panel...\n"
  )

  boxplot_list <- lapply(
    FIGURE2_BOXPLOT_GENES,
    function(
        gene
    ) {

      make_four_group_boxplot(
        gene,
        epithelial_plot_data
      )
    }
  )

  combined_boxplots <- ggpubr::ggarrange(
    plotlist = boxplot_list,
    ncol = 3,
    nrow = 6,
    align = "hv"
  )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "Figure2_CDE_four_group_boxplots.png"
    ),
    plot = combined_boxplots,
    width = 12,
    height = 19,
    units = "in",
    dpi = 400
  )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      "Figure2_CDE_four_group_boxplots.pdf"
    ),
    plot = combined_boxplots,
    width = 12,
    height = 19,
    units = "in"
  )
}


###############################################################################
# 15. Save metadata and session information
###############################################################################

write.csv(
  epithelial_metadata,
  file.path(
    OUTPUT_DIR,
    "Figure2_epithelial_ROI_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)


###############################################################################
# 16. Final console summary
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "GeoMx Figure 2 / early-risk workflow completed\n"
)

cat(
  "============================================================\n"
)

cat(
  "Input: ",
  INPUT_RDS,
  "\n",
  sep = ""
)

cat(
  "Expression scale: log2(q_norm)\n"
)

cat(
  "Epithelial groups: "
)

print(
  table(
    epithelial_metadata$Group
  )
)

cat(
  "Final panel genes: ",
  length(
    FINAL_30_GENES
  ),
  "\n",
  sep = ""
)

cat(
  "Table S2-2 rows: ",
  nrow(
    table_s22
  ),
  "\n",
  sep = ""
)

cat(
  "Four-group boxplots: NO- / NO+ / TO- / TO+\n"
)

cat(
  "Output: ",
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
