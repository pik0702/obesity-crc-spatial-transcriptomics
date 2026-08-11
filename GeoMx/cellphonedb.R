#!/usr/bin/env Rscript

###############################################################################
# GeoMx CellPhoneDB - Figure 1H
#
# Input:
#   epithelial_CD3_means.txt
#   epithelial_CD3_pvalues.txt
#   epithelial_CD3_interaction_scores.txt
#   CD3_epithelial_means.txt
#   CD3_epithelial_pvalues.txt
#   CD3_epithelial_interaction_scores.txt
#   ROI_information_2.csv
#
# Analysis:
#   Ligand-Receptor interactions
#   PValue == 0
#   MeanValue -> dot color
#   Interaction score -> dot size
###############################################################################


###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork"
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
    )
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


###############################################################################
# 2. Paths
###############################################################################

INPUT_DIR <- Sys.getenv(
  "GEOMX_CPDB_PROCESSED_DIR",
  unset = "C:/Users/user/Documents/cellphoneDB/output2"
)

OUTPUT_DIR <- Sys.getenv(
  "GEOMX_CPDB_FIGURE_OUTPUT_DIR",
  unset = "C:/Users/user/Documents/results/GeoMx/06_cellphonedb_figure1H"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 3. Input files
###############################################################################

required_files <- c(
  "epithelial_CD3_means.txt",
  "epithelial_CD3_pvalues.txt",
  "epithelial_CD3_interaction_scores.txt",
  "CD3_epithelial_means.txt",
  "CD3_epithelial_pvalues.txt",
  "CD3_epithelial_interaction_scores.txt",
  "ROI_information_2.csv"
)

missing_files <- required_files[
  !file.exists(
    file.path(
      INPUT_DIR,
      required_files
    )
  )
]

if (
  length(
    missing_files
  ) > 0
) {
  stop(
    "Missing input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


###############################################################################
# 4. Helper functions
###############################################################################

read_cpdb_table <- function(
  filename
) {

  x <- read.delim(
    file.path(
      INPUT_DIR,
      filename
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  valid_names <- !is.na(
    colnames(x)
  ) &
    nzchar(
      colnames(x)
    )

  x[
    ,
    valid_names,
    drop = FALSE
  ]
}


get_patient_columns <- function(
  x
) {

  grep(
    "^Patient_[0-9]+$",
    colnames(x),
    value = TRUE
  )
}


prepare_direction <- function(
  means_file,
  pvalues_file,
  scores_file,
  direction_name
) {

  means_data <- read_cpdb_table(
    means_file
  )

  pvalues_data <- read_cpdb_table(
    pvalues_file
  )

  scores_data <- read_cpdb_table(
    scores_file
  )

  required_columns <- c(
    "interacting_pair",
    "directionality"
  )

  for (
    column_name in
      required_columns
  ) {

    if (
      !column_name %in%
        colnames(
          means_data
        ) ||
        !column_name %in%
          colnames(
            pvalues_data
          ) ||
        !column_name %in%
          colnames(
            scores_data
          )
    ) {
      stop(
        "Required column not found: ",
        column_name
      )
    }
  }


  patient_columns <- get_patient_columns(
    means_data
  )

  if (
    length(
      patient_columns
    ) != 69
  ) {
    warning(
      "Expected 69 patient columns but found ",
      length(
        patient_columns
      ),
      " in ",
      means_file
    )
  }


  if (
    !identical(
      patient_columns,
      get_patient_columns(
        pvalues_data
      )
    ) ||
      !identical(
        patient_columns,
        get_patient_columns(
          scores_data
        )
      )
  ) {
    stop(
      "Patient columns differ among means, p-values, and interaction scores."
    )
  }


  if (
    !identical(
      as.character(
        means_data$interacting_pair
      ),
      as.character(
        pvalues_data$interacting_pair
      )
    ) ||
      !identical(
        as.character(
          means_data$interacting_pair
        ),
        as.character(
          scores_data$interacting_pair
        )
      )
  ) {
    stop(
      "Interaction order differs among means, p-values, and scores."
    )
  }


  ligand_receptor_rows <- means_data$directionality ==
    "Ligand-Receptor"


  means_long <- means_data[
    ligand_receptor_rows,
    c(
      "interacting_pair",
      patient_columns
    ),
    drop = FALSE
  ] |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patient_columns
      ),
      names_to = "Patient",
      values_to = "MeanValue"
    )


  pvalues_long <- pvalues_data[
    ligand_receptor_rows,
    c(
      "interacting_pair",
      patient_columns
    ),
    drop = FALSE
  ] |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patient_columns
      ),
      names_to = "Patient",
      values_to = "PValue"
    )


  scores_long <- scores_data[
    ligand_receptor_rows,
    c(
      "interacting_pair",
      patient_columns
    ),
    drop = FALSE
  ] |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patient_columns
      ),
      names_to = "Patient",
      values_to = "Score"
    )


  out <- means_long |>
    dplyr::left_join(
      pvalues_long,
      by = c(
        "interacting_pair",
        "Patient"
      )
    ) |>
    dplyr::left_join(
      scores_long,
      by = c(
        "interacting_pair",
        "Patient"
      )
    ) |>
    dplyr::mutate(
      MeanValue = as.numeric(
        MeanValue
      ),
      PValue = as.numeric(
        PValue
      ),
      Score = as.numeric(
        Score
      ),
      PatientID = sub(
        "^Patient_",
        "",
        Patient
      ),
      Direction = direction_name
    )

  message(
    direction_name,
    ": ",
    dplyr::n_distinct(
      out$interacting_pair
    ),
    " Ligand-Receptor interactions; ",
    dplyr::n_distinct(
      out$PatientID
    ),
    " patients"
  )

  out
}


###############################################################################
# 5. Read CellPhoneDB matrices
###############################################################################

tumor_to_immune <- prepare_direction(
  means_file =
    "epithelial_CD3_means.txt",
  pvalues_file =
    "epithelial_CD3_pvalues.txt",
  scores_file =
    "epithelial_CD3_interaction_scores.txt",
  direction_name =
    "Tumor - Immune"
)

immune_to_tumor <- prepare_direction(
  means_file =
    "CD3_epithelial_means.txt",
  pvalues_file =
    "CD3_epithelial_pvalues.txt",
  scores_file =
    "CD3_epithelial_interaction_scores.txt",
  direction_name =
    "Immune - Tumor"
)

cpdb_data <- dplyr::bind_rows(
  tumor_to_immune,
  immune_to_tumor
)


###############################################################################
# 6. Patient annotation
###############################################################################

roi_info <- read.csv(
  file.path(
    INPUT_DIR,
    "ROI_information_2.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

valid_names <- !is.na(
  colnames(
    roi_info
  )
) &
  nzchar(
    colnames(
      roi_info
    )
  )

roi_info <- roi_info[
  ,
  valid_names,
  drop = FALSE
]

required_annotation_columns <- c(
  "PatientID",
  "type",
  "class"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(
    roi_info
  )
)

if (
  length(
    missing_annotation_columns
  ) > 0
) {
  stop(
    "Missing ROI annotation column(s): ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
}


patient_annotation <- roi_info |>
  dplyr::filter(
    class ==
      "Tumor"
  ) |>
  dplyr::transmute(
    PatientID = as.character(
      PatientID
    ),
    Obesity = dplyr::case_when(
      tolower(
        trimws(
          as.character(
            type
          )
        )
      ) == "obese" ~
        "Obese",

      TRUE ~
        "Non-Obese"
    )
  ) |>
  dplyr::distinct()


duplicate_annotation <- patient_annotation |>
  dplyr::count(
    PatientID,
    name = "N"
  ) |>
  dplyr::filter(
    N > 1
  )

if (
  nrow(
    duplicate_annotation
  ) > 0
) {
  stop(
    "Conflicting obesity annotation detected for patient(s): ",
    paste(
      duplicate_annotation$PatientID,
      collapse = ", "
    )
  )
}


cpdb_data <- cpdb_data |>
  dplyr::left_join(
    patient_annotation,
    by = "PatientID"
  )

if (
  anyNA(
    cpdb_data$Obesity
  )
) {

  missing_ids <- unique(
    cpdb_data$PatientID[
      is.na(
        cpdb_data$Obesity
      )
    ]
  )

  stop(
    "Missing obesity annotation for patient(s): ",
    paste(
      missing_ids,
      collapse = ", "
    )
  )
}


patient_summary <- cpdb_data |>
  dplyr::distinct(
    PatientID,
    Obesity
  ) |>
  dplyr::count(
    Obesity,
    name = "N_patients"
  )

print(
  patient_summary
)


###############################################################################
# 7. Figure 1H interactions
###############################################################################

figure_interactions <- data.frame(
  Direction = c(
    rep(
      "Tumor - Immune",
      5
    ),
    rep(
      "Immune - Tumor",
      10
    )
  ),

  interacting_pair = c(
    "PVR_TIGIT",
    "AREG_EGFR",
    "LTB_LTBR",
    "IGF1_IGF1R",
    "EFNA4_EPHA1",
    "WNT8A_FZD5_LRP6",
    "MICA_NKG2D_II_receptor",
    "JAG1_NOTCH3",
    "IL25_IL17_receptor_AB",
    "GAS6_AXL",
    "WNT10A_FZD5_LRP6",
    "THBS2_CD36",
    "IGF2_IGF1R",
    "HLA-F_LILRB2",
    "CTSG_FPR1"
  ),

  Display = c(
    "PVR - TIGIT",
    "AREG - EGFR",
    "LTB - LTBR",
    "IGF1 - IGF1R",
    "EFNA4 - EPHA1",
    "WNT8A - FZD5 - LRP6",
    "MICA - NKG2D II receptor",
    "JAG1 - NOTCH3",
    "IL25 - IL17 receptor AB",
    "GAS6 - AXL",
    "WNT10A - FZD5 - LRP6",
    "THBS2 - CD36",
    "IGF2 - IGF1R",
    "HLA-F - LILRB2",
    "CTSG - FPR1"
  ),

  stringsAsFactors = FALSE
)


missing_interactions <- figure_interactions |>
  dplyr::anti_join(
    cpdb_data |>
      dplyr::distinct(
        Direction,
        interacting_pair
      ),
    by = c(
      "Direction",
      "interacting_pair"
    )
  )

if (
  nrow(
    missing_interactions
  ) > 0
) {
  stop(
    "Figure 1H interaction(s) missing:\n",
    paste(
      missing_interactions$interacting_pair,
      collapse = "\n"
    )
  )
}


figure1h_data <- cpdb_data |>
  dplyr::inner_join(
    figure_interactions,
    by = c(
      "Direction",
      "interacting_pair"
    )
  ) |>
  dplyr::filter(
    PValue ==
      0
  )


###############################################################################
# 8. Patient order
###############################################################################

patient_order <- cpdb_data |>
  dplyr::distinct(
    PatientID,
    Obesity
  ) |>
  dplyr::mutate(
    PatientNumber = suppressWarnings(
      as.numeric(
        PatientID
      )
    ),
    ObesityOrder = ifelse(
      Obesity ==
        "Non-Obese",
      1,
      2
    )
  ) |>
  dplyr::arrange(
    ObesityOrder,
    PatientNumber
  ) |>
  dplyr::pull(
    PatientID
  )


figure1h_data <- figure1h_data |>
  dplyr::mutate(
    PatientID = factor(
      PatientID,
      levels =
        patient_order
    ),
    Obesity = factor(
      Obesity,
      levels = c(
        "Non-Obese",
        "Obese"
      )
    )
  )


###############################################################################
# 9. Export plotting data
###############################################################################

write.csv(
  figure1h_data,
  file.path(
    OUTPUT_DIR,
    "Figure1H_plot_data.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  patient_annotation |>
    dplyr::filter(
      PatientID %in%
        patient_order
    ) |>
    dplyr::arrange(
      match(
        PatientID,
        patient_order
      )
    ),
  file.path(
    OUTPUT_DIR,
    "Figure1H_patient_annotation.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  patient_summary,
  file.path(
    OUTPUT_DIR,
    "Figure1H_patient_counts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


bubble_summary <- figure1h_data |>
  dplyr::count(
    Direction,
    interacting_pair,
    Obesity,
    name = "N_bubbles"
  )

write.csv(
  bubble_summary,
  file.path(
    OUTPUT_DIR,
    "Figure1H_bubble_counts.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 10. Plot
###############################################################################

tumor_immune_levels <- figure_interactions |>
  dplyr::filter(
    Direction ==
      "Tumor - Immune"
  ) |>
  dplyr::pull(
    Display
  )

immune_tumor_levels <- figure_interactions |>
  dplyr::filter(
    Direction ==
      "Immune - Tumor"
  ) |>
  dplyr::pull(
    Display
  )


make_panel <- function(
  data,
  direction_name,
  display_levels
) {

  plot_data <- data |>
    dplyr::filter(
      Direction ==
        direction_name
    ) |>
    dplyr::mutate(
      Display = factor(
        Display,
        levels = rev(
          display_levels
        )
      )
    )

  ggplot(
    plot_data,
    aes(
      x = PatientID,
      y = Display
    )
  ) +
    geom_point(
      aes(
        size = Score,
        color = MeanValue
      )
    ) +
    facet_grid(
      . ~ Obesity,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_color_gradient(
      name = "Mean",
      low = "lightblue",
      high = "red"
    ) +
    scale_size_continuous(
      limits = c(
        0,
        100
      ),
      breaks = c(
        0,
        25,
        50,
        75,
        100
      ),
      range = c(
        0.5,
        5
      ),
      name = "Interaction\nScore"
    ) +
    labs(
      title = direction_name,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(
      base_size = 11
    ) +
    theme(
      axis.text.x =
        element_blank(),
      axis.ticks.x =
        element_blank(),
      axis.text.y =
        element_text(
          color = "black",
          size = 10
        ),
      panel.grid.minor =
        element_blank(),
      strip.background =
        element_rect(
          fill = "white",
          color = "black"
        ),
      strip.text =
        element_text(
          color = "black",
          size = 10
        ),
      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold",
          size = 13
        ),
      legend.position =
        "right"
    )
}


plot_tumor_immune <- make_panel(
  figure1h_data,
  "Tumor - Immune",
  tumor_immune_levels
)

plot_immune_tumor <- make_panel(
  figure1h_data,
  "Immune - Tumor",
  immune_tumor_levels
)


figure1h <- plot_tumor_immune /
  plot_immune_tumor +
  patchwork::plot_layout(
    heights = c(
      1,
      1.6
    )
  )


ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure1H_CellPhoneDB.pdf"
  ),
  plot = figure1h,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure1H_CellPhoneDB.png"
  ),
  plot = figure1h,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)


###############################################################################
# 11. Summary
###############################################################################

cat(
  "\nFigure 1H data rows:",
  nrow(
    figure1h_data
  ),
  "\n"
)

cat(
  "Patients:",
  dplyr::n_distinct(
    figure1h_data$PatientID
  ),
  "\n"
)

cat(
  "Output:",
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  ),
  "\n"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTPUT_DIR,
    "sessionInfo_R.txt"
  )
)
