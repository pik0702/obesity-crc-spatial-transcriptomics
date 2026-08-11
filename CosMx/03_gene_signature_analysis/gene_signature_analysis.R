#!/usr/bin/env Rscript

###############################################################################
# 03_CosMx_gene_signature_analysis.R
#
# CosMx gene-level signature analysis
# Obesity-associated colorectal cancer study
#
# Final annotation column
# -----------------------
# Major cell type : cell_type_major
#
# This script intentionally does NOT calculate composite module scores.
# AddModuleScore or other per-cell scoring methods are not used.
#
# Workflow
# --------
# 1. Load final annotated CosMx Seurat object
# 2. Subset Tumor Cells and Myeloids using cell_type_major
# 3. Evaluate GeoMx-derived biomarker genes in Tumor Cells and Myeloids
# 4. Compare gene expression between Obese and Non-obese groups
# 5. Evaluate M1/M2 polarization marker genes in Myeloids
# 6. Save gene-level expression summaries and statistical results
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- "path/to/cosmx_final.rds"

OUTPUT_DIR <- "results/CosMx/03_gene_signature_analysis"

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    OUTPUT_DIR,
    "plots"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. Packages
###############################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


###############################################################################
# 2. Load final annotated CosMx object
###############################################################################

cosmx <- readRDS(
  INPUT_RDS
)


###############################################################################
# 3. Final metadata columns
###############################################################################

MAJOR_COLUMN <- "cell_type_major"

OBESITY_COLUMN <- "obese_status"


###############################################################################
# 4. GeoMx-derived biomarker gene panel
#
# Gene panel used in the original CosMx validation analysis.
###############################################################################

geomx_biomarker_genes <- c(
  "NES",
  "CMSS1",
  "BGN",
  "CCL24",
  "PTPN13",
  "FGFRL1",
  "MYL9",
  "RAB22A",
  "MCL1",
  "AGR2",
  "COX5A",
  "SPINK1",
  "TMPRSS2",
  "MUC2"
)


###############################################################################
# 5. Tumor Cells and Myeloids
###############################################################################

tumor_cells <- subset(
  cosmx,
  subset = cell_type_major == "Tumor Cells"
)

myeloids <- subset(
  cosmx,
  subset = cell_type_major == "Myeloids"
)


###############################################################################
# 6. Keep genes available in the CosMx dataset
###############################################################################

geomx_genes_use <- intersect(
  geomx_biomarker_genes,
  rownames(cosmx)
)


###############################################################################
# 7. Helper function: gene-level Obese vs Non-obese analysis
#
# No composite score is calculated.
# Each gene is analyzed independently.
###############################################################################

analyze_gene_panel <- function(
    object,
    genes,
    cell_group
) {

  genes_use <- intersect(
    genes,
    rownames(object)
  )

  expression_data <- FetchData(
    object,
    vars = c(
      genes_use,
      OBESITY_COLUMN
    ),
    layer = "data"
  )

  expression_data$cell <- rownames(
    expression_data
  )

  expression_long <- expression_data %>%
    pivot_longer(
      cols = all_of(
        genes_use
      ),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    mutate(
      cell_group = cell_group
    )


  ###########################################################################
  # Mean and median expression by obesity status
  ###########################################################################

  expression_summary <- expression_long %>%
    group_by(
      cell_group,
      gene,
      .data[[OBESITY_COLUMN]]
    ) %>%
    summarise(
      mean_expression = mean(
        expression,
        na.rm = TRUE
      ),
      median_expression = median(
        expression,
        na.rm = TRUE
      ),
      n_cells = n(),
      .groups = "drop"
    )


  ###########################################################################
  # Gene-wise Wilcoxon test
  ###########################################################################

  test_results <- expression_long %>%
    group_by(
      cell_group,
      gene
    ) %>%
    group_modify(
      ~ {

        dat <- .x %>%
          filter(
            !is.na(
              .data[[OBESITY_COLUMN]]
            )
          )

        groups_present <- unique(
          dat[[OBESITY_COLUMN]]
        )

        if (
          length(groups_present) < 2
        ) {

          return(
            data.frame(
              p_value = NA_real_
            )
          )
        }

        test_formula <- reformulate(
          OBESITY_COLUMN,
          response = "expression"
        )

        data.frame(
          p_value = wilcox.test(
            test_formula,
            data = dat,
            exact = FALSE
          )$p.value
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      p_adjust = p.adjust(
        p_value,
        method = "BH"
      )
    )


  list(
    expression = expression_long,
    summary = expression_summary,
    statistics = test_results
  )
}


###############################################################################
# 8. GeoMx-derived biomarker analysis in Tumor Cells
###############################################################################

geomx_tumor <- analyze_gene_panel(
  object = tumor_cells,
  genes = geomx_genes_use,
  cell_group = "Tumor Cells"
)


###############################################################################
# 9. GeoMx-derived biomarker analysis in Myeloids
###############################################################################

geomx_myeloid <- analyze_gene_panel(
  object = myeloids,
  genes = geomx_genes_use,
  cell_group = "Myeloids"
)


###############################################################################
# 10. Combine GeoMx biomarker results
###############################################################################

geomx_expression_summary <- bind_rows(
  geomx_tumor$summary,
  geomx_myeloid$summary
)

geomx_statistics <- bind_rows(
  geomx_tumor$statistics,
  geomx_myeloid$statistics
)


write.csv(
  geomx_expression_summary,
  file.path(
    OUTPUT_DIR,
    "GeoMx_biomarker_expression_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  geomx_statistics,
  file.path(
    OUTPUT_DIR,
    "GeoMx_biomarker_obese_vs_nonobese_statistics.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 11. GeoMx-derived biomarker DotPlot
#
# Tumor Cells and Myeloids are shown separately.
###############################################################################

tumor_cells$plot_group <- paste0(
  tumor_cells[[OBESITY_COLUMN]][, 1],
  "_Tumor"
)

myeloids$plot_group <- paste0(
  myeloids[[OBESITY_COLUMN]][, 1],
  "_Myeloids"
)


combined_for_plot <- merge(
  tumor_cells,
  y = myeloids
)


p_geomx <- DotPlot(
  combined_for_plot,
  features = geomx_genes_use,
  group.by = "plot_group"
) +
  RotatedAxis() +
  labs(
    x = "",
    y = "",
    title = "GeoMx-derived biomarker genes"
  ) +
  theme_minimal()


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "GeoMx_biomarker_Tumor_Myeloid_DotPlot.jpeg"
  ),
  p_geomx,
  width = 11,
  height = 5,
  dpi = 300
)


###############################################################################
# 12. M1/M2 polarization marker panel
#
# Gene-level polarization markers recovered from the original CosMx analysis.
# These are analyzed individually; no M1 or M2 composite score is calculated.
###############################################################################

M1_markers <- c(
  "CD80",
  "CD86",
  "TNF",
  "IL1B",
  "CXCL9",
  "CXCL10"
)

M2_markers <- c(
  "CD163",
  "ARG1",
  "MRC1",
  "IL10",
  "TGFB1",
  "CCL18"
)

polarization_genes <- c(
  M1_markers,
  M2_markers
)


###############################################################################
# 13. M1/M2 polarization gene-level analysis in Myeloids
###############################################################################

polarization_result <- analyze_gene_panel(
  object = myeloids,
  genes = polarization_genes,
  cell_group = "Myeloids"
)


polarization_summary <- polarization_result$summary %>%
  mutate(
    polarization = case_when(
      gene %in% M1_markers ~ "M1",
      gene %in% M2_markers ~ "M2",
      TRUE ~ NA_character_
    )
  )


polarization_statistics <- polarization_result$statistics %>%
  mutate(
    polarization = case_when(
      gene %in% M1_markers ~ "M1",
      gene %in% M2_markers ~ "M2",
      TRUE ~ NA_character_
    )
  )


write.csv(
  polarization_summary,
  file.path(
    OUTPUT_DIR,
    "M1_M2_polarization_gene_expression_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  polarization_statistics,
  file.path(
    OUTPUT_DIR,
    "M1_M2_polarization_obese_vs_nonobese_statistics.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 14. M1/M2 polarization DotPlot
###############################################################################

myeloid_genes_use <- intersect(
  polarization_genes,
  rownames(myeloids)
)


p_polarization <- DotPlot(
  myeloids,
  features = list(
    M1 = intersect(
      M1_markers,
      myeloid_genes_use
    ),
    M2 = intersect(
      M2_markers,
      myeloid_genes_use
    )
  ),
  group.by = OBESITY_COLUMN
) +
  RotatedAxis() +
  labs(
    x = "",
    y = "",
    title = "M1/M2 polarization markers"
  ) +
  theme_minimal()


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "M1_M2_polarization_DotPlot.jpeg"
  ),
  p_polarization,
  width = 9,
  height = 4,
  dpi = 300
)
