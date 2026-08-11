#!/usr/bin/env Rscript

###############################################################################
# 04_CosMx_spatial_niche_analysis.R
#
# CosMx spatial niche analysis
# Obesity-associated colorectal cancer study
#
# Three spatial compartments were considered:
#   - Intratumoral
#   - Intertumoral
#   - Adipose-associated
#
# The code below includes one representative crop from the original workflow.
# The same procedure was applied to the remaining manually selected regions.
#
# Workflow
# --------
# 1. Load annotated CosMx Seurat object
# 2. Define manually selected crop regions
# 3. Crop regions and calculate area (mm^2)
# 4. Quantify major cell-type density (cells/mm^2)
# 5. Quantify SPP1+ and Pro-inflammatory macrophage density
# 6. Calculate SPP1+ / Pro-inflammatory ratio by spatial compartment
# 7. Extract cell coordinates
# 8. Calculate SPP1+ macrophage nearest-neighbor distance to tumor cells
# 9. Compare obese vs non-obese samples
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- "path/to/CosMx_annotated_object.rds"

OUTPUT_DIR <- "results/CosMx/04_spatial_niche_analysis"

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
# 1. Packages
###############################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FNN)


###############################################################################
# 2. Load annotated CosMx Seurat object
###############################################################################

cosmx <- readRDS(INPUT_RDS)


###############################################################################
# 3. Metadata columns used in the original Figure 5 workflow
###############################################################################

# Detailed annotation used for SPP1+ and Pro-inflammatory macrophages
SUBTYPE_COLUMN <- "cell_type_figure5_032825"

# Annotation used for major cell-type grouping
MAJOR_COLUMN <- "cell_type_figure5_combine"

OBESITY_COLUMN <- "obese_status"
PATIENT_COLUMN <- "patient"


###############################################################################
# 4. Manual spatial-region definition
#
# Regions were selected by inspection of tissue pathology together with
# local enrichment of SPP1+ and Pro-inflammatory macrophages.
#
# Representative example from the original workflow:
#   s1p1_t1
#   image: colon.TMA_JHK_1
#   x = 14.20 - 14.57 mm
#   y =  1.55 -  1.68 mm
#
# Add the remaining manually selected regions as additional rows using
# the same format. Set spatial_niche to one of:
#   "Intratumoral"
#   "Intertumoral"
#   "Adipose-associated"
###############################################################################

selected_regions <- data.frame(
  spatial_field = c(
    "s1p1_t1"
  ),
  image = c(
    "colon.TMA_JHK_1"
  ),
  spatial_niche = c(
    "Intratumoral"
  ),
  x_min = c(
    14.20
  ),
  x_max = c(
    14.57
  ),
  y_min = c(
    1.55
  ),
  y_max = c(
    1.68
  ),
  stringsAsFactors = FALSE
)


###############################################################################
# 5. Crop selected regions and calculate region area
###############################################################################

selected_regions <- selected_regions %>%
  mutate(
    Area_mm2 = abs(x_max - x_min) *
      abs(y_max - y_min)
  )


cosmx@meta.data$spatial_field <- NA_character_
cosmx@meta.data$spatial_niche <- NA_character_


for (i in seq_len(nrow(selected_regions))) {

  region <- selected_regions[i, ]

  crop_object <- Crop(
    cosmx[[region$image]],
    x = c(
      region$x_min,
      region$x_max
    ),
    y = c(
      region$y_min,
      region$y_max
    )
  )

  cosmx[[region$spatial_field]] <- crop_object

  DefaultBoundary(
    cosmx[[region$spatial_field]]
  ) <- "segmentation"

  cropped_cells <- Cells(
    cosmx[[region$spatial_field]]
  )

  cosmx@meta.data$spatial_field[
    rownames(cosmx@meta.data) %in% cropped_cells
  ] <- region$spatial_field

  cosmx@meta.data$spatial_niche[
    rownames(cosmx@meta.data) %in% cropped_cells
  ] <- region$spatial_niche
}


write.csv(
  selected_regions,
  file.path(
    OUTPUT_DIR,
    "manual_spatial_regions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 6. Major cell-type grouping
#
# This follows the grouping used in the original Figure 5 workflow.
###############################################################################

meta <- cosmx@meta.data %>%
  tibble::rownames_to_column("cell")


meta$major_cell_type <- meta[[MAJOR_COLUMN]]


# Epithelial/Stromal labels were grouped with Tumor Cells
meta$major_cell_type[
  meta$major_cell_type %in% c(
    "Epithelial Cells",
    "Stromal Cells"
  )
] <- "Tumor Cells"


# Collapse detailed immune annotations into major cell types
meta <- meta %>%
  mutate(
    major_cell_type = case_when(

      major_cell_type %in% c(
        "CD19+CD20+ B",
        "IgG+ Plasma",
        "IgA+ Plasma"
      ) ~ "B cells",

      major_cell_type %in% c(
        "CD4+ T cells",
        "CD8+ T cells",
        "NK cells",
        "Regulatory T cells"
      ) ~ "T cells",

      major_cell_type %in% c(
        "SPP1+",
        "Pro-inflammatory",
        "Myloids"
      ) ~ "Myeloids",

      TRUE ~ major_cell_type
    )
  )


###############################################################################
# 7. Major cell-type density in manually selected regions
#
# Density = number of cells / cropped area (mm^2)
###############################################################################

major_density <- meta %>%
  filter(
    !is.na(spatial_field),
    !is.na(.data[[OBESITY_COLUMN]])
  ) %>%
  count(
    spatial_field,
    spatial_niche,
    !!sym(OBESITY_COLUMN),
    major_cell_type,
    name = "count"
  ) %>%
  left_join(
    selected_regions %>%
      select(
        spatial_field,
        Area_mm2
      ),
    by = "spatial_field"
  ) %>%
  mutate(
    Count_per_mm2 = count / Area_mm2
  )


write.csv(
  major_density,
  file.path(
    OUTPUT_DIR,
    "major_celltype_density_per_mm2.csv"
  ),
  row.names = FALSE
)


# Obese vs non-obese density comparison
p_major <- ggplot(
  major_density,
  aes(
    x = .data[[OBESITY_COLUMN]],
    y = Count_per_mm2,
    fill = .data[[OBESITY_COLUMN]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ major_cell_type,
    scales = "free_y"
  ) +
  labs(
    x = "Obesity Status",
    y = "Cell Density (Count per mm²)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "major_celltype_density_obese_vs_nonobese.jpeg"
  ),
  plot = p_major,
  width = 12,
  height = 8,
  dpi = 300
)


###############################################################################
# 8. SPP1+ and Pro-inflammatory macrophage density
###############################################################################

macrophage_density <- meta %>%
  filter(
    !is.na(spatial_field),
    .data[[SUBTYPE_COLUMN]] %in% c(
      "SPP1+",
      "Pro-inflammatory"
    ),
    !is.na(.data[[OBESITY_COLUMN]])
  ) %>%
  count(
    spatial_field,
    spatial_niche,
    !!sym(OBESITY_COLUMN),
    !!sym(SUBTYPE_COLUMN),
    name = "count"
  ) %>%
  left_join(
    selected_regions %>%
      select(
        spatial_field,
        Area_mm2
      ),
    by = "spatial_field"
  ) %>%
  mutate(
    Count_per_mm2 = count / Area_mm2
  )


write.csv(
  macrophage_density,
  file.path(
    OUTPUT_DIR,
    "SPP1_ProInflammatory_density_per_mm2.csv"
  ),
  row.names = FALSE
)


p_macrophage <- ggplot(
  macrophage_density,
  aes(
    x = .data[[OBESITY_COLUMN]],
    y = Count_per_mm2,
    fill = .data[[OBESITY_COLUMN]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  facet_grid(
    spatial_niche ~ .data[[SUBTYPE_COLUMN]],
    scales = "free_y"
  ) +
  labs(
    x = "Obesity Status",
    y = "Cell Density (Count per mm²)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "SPP1_ProInflammatory_density_obese_vs_nonobese.jpeg"
  ),
  plot = p_macrophage,
  width = 8,
  height = 8,
  dpi = 300
)


###############################################################################
# 9. SPP1+ / Pro-inflammatory ratio
#
# Ratio is calculated from density (cells/mm^2), following the original
# workflow. Because both cell types are measured in the same cropped region,
# this is equivalent to the corresponding cell-count ratio.
#
# Ratios are evaluated separately for:
#   - Intratumoral
#   - Intertumoral
#   - Adipose-associated
###############################################################################

ratio_data <- macrophage_density %>%
  select(
    spatial_field,
    spatial_niche,
    !!sym(OBESITY_COLUMN),
    !!sym(SUBTYPE_COLUMN),
    Count_per_mm2
  ) %>%
  pivot_wider(
    names_from = !!sym(SUBTYPE_COLUMN),
    values_from = Count_per_mm2,
    values_fill = 0
  ) %>%
  mutate(
    SPP1_to_Pro_ratio = ifelse(
      `Pro-inflammatory` > 0,
      `SPP1+` / `Pro-inflammatory`,
      NA_real_
    )
  )


write.csv(
  ratio_data,
  file.path(
    OUTPUT_DIR,
    "SPP1_to_ProInflammatory_ratio_by_niche.csv"
  ),
  row.names = FALSE
)


p_ratio <- ggplot(
  ratio_data,
  aes(
    x = .data[[OBESITY_COLUMN]],
    y = SPP1_to_Pro_ratio,
    fill = .data[[OBESITY_COLUMN]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ spatial_niche,
    scales = "free_y"
  ) +
  labs(
    x = "Obesity Status",
    y = "SPP1+ / Pro-inflammatory ratio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "SPP1_to_ProInflammatory_ratio_by_niche.jpeg"
  ),
  plot = p_ratio,
  width = 9,
  height = 4,
  dpi = 300
)


###############################################################################
# 10. Extract cell coordinates
#
# Original CosMx slides used for the coordinate-based analysis.
###############################################################################

slide_images <- c(
  "colon.TMA_JHK_1",
  "colon.TMA_JHK_2"
)


coordinate_list <- lapply(
  slide_images,
  function(image_name) {

    coords <- GetTissueCoordinates(
      cosmx,
      image = image_name,
      type = "centroids"
    )

    coords <- as.data.frame(coords)

    coords$cell <- as.character(
      coords$cell
    )

    coords$image <- image_name

    coords
  }
)


coords_all <- bind_rows(
  coordinate_list
)


meta_coordinates <- cosmx@meta.data %>%
  tibble::rownames_to_column("cell")


coords_full <- left_join(
  coords_all,
  meta_coordinates,
  by = "cell"
)


###############################################################################
# 11. SPP1+ macrophage nearest-neighbor distance to tumor cells
#
# For each patient:
#   query = SPP1+ macrophage coordinates
#   data  = Tumor-cell coordinates
#   k     = 1
#
# Only SPP1+ macrophages are retained in this public workflow.
###############################################################################

patients <- unique(
  coords_full[[PATIENT_COLUMN]]
)


all_spp1_distances <- list()


for (pt in patients) {

  data_pt <- coords_full %>%
    filter(
      .data[[PATIENT_COLUMN]] == pt
    )


  tumor_data <- data_pt %>%
    filter(
      .data[[SUBTYPE_COLUMN]] == "Tumor Cells"
    )


  spp1_data <- data_pt %>%
    filter(
      .data[[SUBTYPE_COLUMN]] == "SPP1+"
    )


  if (
    nrow(tumor_data) < 1 ||
    nrow(spp1_data) < 1
  ) {
    next
  }


  tumor_coords <- as.matrix(
    tumor_data %>%
      select(x, y)
  )


  spp1_coords <- as.matrix(
    spp1_data %>%
      select(x, y)
  )


  nn_result <- FNN::get.knnx(
    data = tumor_coords,
    query = spp1_coords,
    k = 1
  )


  spp1_data$nn_distance <- nn_result$nn.dist[, 1]


  all_spp1_distances[[
    length(all_spp1_distances) + 1
  ]] <- spp1_data %>%
    select(
      cell,
      all_of(PATIENT_COLUMN),
      all_of(OBESITY_COLUMN),
      nn_distance
    )
}


spp1_nn <- bind_rows(
  all_spp1_distances
)


write.csv(
  spp1_nn,
  file.path(
    OUTPUT_DIR,
    "SPP1_nearest_tumor_distance.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 12. Obese vs non-obese comparison of SPP1+ NN distance
###############################################################################

p_nn <- ggplot(
  spp1_nn,
  aes(
    x = .data[[OBESITY_COLUMN]],
    y = nn_distance,
    fill = .data[[OBESITY_COLUMN]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,
    alpha = 0.4
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  labs(
    x = "Obesity Status",
    y = "Nearest Tumor Cell Distance (mm)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "SPP1_nearest_tumor_distance_obese_vs_nonobese.jpeg"
  ),
  plot = p_nn,
  width = 5,
  height = 5,
  dpi = 300
)


###############################################################################
# 13. Save object with spatial-region metadata
###############################################################################

saveRDS(
  cosmx,
  file.path(
    OUTPUT_DIR,
    "CosMx_spatial_niche_annotated.rds"
  )
)
