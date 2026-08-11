#!/usr/bin/env Rscript

###############################################################################
# GeoMx DSP -> CellPhoneDB -> Figure 1H
#
# Steps:
#   1. Load the Q3-normalized GeoMx object.
#   2. Select tumor CD3, Epithelial, and Adipose ROIs.
#   3. Transform Q3 expression as log2(Q3).
#   4. Write patient-level CellPhoneDB metadata and expression inputs.
#   5. Run CellPhoneDB statistical analysis for each patient.
#   6. Merge Epithelial -> CD3 and CD3 -> Epithelial results.
#   7. Filter ligand-receptor interactions with PValue == 0.
#   8. Draw Figure 1H.
###############################################################################


###############################################################################
# 0. Configuration
###############################################################################

PREPROCESSED_RDS <- Sys.getenv(
  "GEOMX_PREPROCESSED_RDS",
  unset = "C:/Users/user/Documents/results/GeoMx/01_preprocessing_QC_normalization/GeoMx_preprocessed_Q3_normalized.rds"
)

OUTPUT_DIR <- Sys.getenv(
  "GEOMX_CPDB_OUTPUT_DIR",
  unset = "C:/Users/user/Documents/results/GeoMx/05_cellphonedb_full"
)

PYTHON_SCRIPT <- Sys.getenv(
  "GEOMX_CPDB_PYTHON_SCRIPT",
  unset = "C:/Users/user/Documents/geomx/05_run_cellphonedb.py"
)

PYTHON_BIN <- Sys.getenv(
  "GEOMX_CPDB_PYTHON",
  unset = "C:/Users/user/cpdb310/Scripts/python.exe"
)

CPDB_FILE <- Sys.getenv(
  "GEOMX_CPDB_FILE",
  unset = "C:/Users/user/Documents/cellphoneDB/cpdb/cellphonedb.zip"
)


FIGURE1H_INPUT_DIR <- Sys.getenv(
  "GEOMX_CPDB_FIGURE1H_INPUT_DIR",
  unset = "C:/Users/user/Documents/cellphoneDB/output2"
)

RUN_CPDB <- tolower(
  Sys.getenv(
    "GEOMX_CPDB_RUN",
    unset = "true"
  )
) %in% c(
  "true",
  "1",
  "yes",
  "y"
)

PATIENT_SET_MODE <- Sys.getenv(
  "GEOMX_CPDB_PATIENT_SET_MODE",
  unset = "figure1h"
)

# Accept the previous setting name as an alias.
if (
  PATIENT_SET_MODE %in% c(
    "figure1h_legacy",
    "legacy"
  )
) {
  PATIENT_SET_MODE <- "figure1h"
}

CPDB_ITERATIONS <- as.integer(
  Sys.getenv(
    "GEOMX_CPDB_ITERATIONS",
    unset = "1000"
  )
)

CPDB_THRESHOLD <- as.numeric(
  Sys.getenv(
    "GEOMX_CPDB_THRESHOLD",
    unset = "0.1"
  )
)

CPDB_THREADS <- as.integer(
  Sys.getenv(
    "GEOMX_CPDB_THREADS",
    unset = "4"
  )
)

CPDB_COUNTS_DATA <- Sys.getenv(
  "GEOMX_CPDB_COUNTS_DATA",
  unset = "hgnc_symbol"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

INPUT_DIR <- file.path(
  OUTPUT_DIR,
  "inputs"
)

META_DIR <- file.path(
  INPUT_DIR,
  "ROI_table_per_patients"
)

EXPR_DIR <- file.path(
  INPUT_DIR,
  "Expr_Matrix"
)


RAW_DIR <- file.path(
  OUTPUT_DIR,
  "raw_patient_outputs"
)

MERGED_DIR <- file.path(
  OUTPUT_DIR,
  "merged"
)

FIGURE_DIR <- file.path(
  OUTPUT_DIR,
  "figures"
)

for (d in c(
  INPUT_DIR,
  META_DIR,
  EXPR_DIR,
  RAW_DIR,
  MERGED_DIR,
  FIGURE_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


###############################################################################
# 1. Required R packages
###############################################################################

required_packages <- c(
  "Biobase",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales"
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
    "Missing required R package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})


###############################################################################
# 2. Helper functions
###############################################################################

clean_roi_id <- function(x) {

  x <- basename(
    as.character(x)
  )

  x <- sub(
    "\\.dcc$",
    "",
    x,
    ignore.case = TRUE
  )

  x
}


resolve_column <- function(
  data,
  candidates,
  label
) {

  hit <- candidates[
    candidates %in% colnames(data)
  ]

  if (length(hit) == 0) {
    stop(
      "Could not find metadata column for ",
      label,
      ". Tried: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  hit[1]
}


standardize_class <- function(x) {

  x <- trimws(
    as.character(x)
  )

  dplyr::case_when(
    tolower(x) == "normal" ~ "Normal",
    tolower(x) %in% c(
      "tumor",
      "tumour"
    ) ~ "Tumor",
    TRUE ~ x
  )
}


standardize_obesity <- function(x) {

  x <- trimws(
    as.character(x)
  )

  low <- tolower(x)

  dplyr::case_when(
    low %in% c(
      "normal",
      "non-obese",
      "nonobese",
      "non obese",
      "non_obese"
    ) ~ "Normal",

    low == "obese" ~ "Obese",

    TRUE ~ x
  )
}


standardize_roi_type <- function(
  roi_type,
  segment = NULL
) {

  x <- trimws(
    as.character(roi_type)
  )

  if (!is.null(segment)) {

    segment <- trimws(
      as.character(segment)
    )

    x[
      is.na(x) |
        x == ""
    ] <- segment[
      is.na(x) |
        x == ""
    ]
  }

  dplyr::case_when(
    tolower(x) %in% c(
      "cd3",
      "immune",
      "immune cell"
    ) ~ "CD3",

    tolower(x) %in% c(
      "adipose",
      "adipose tissue"
    ) ~ "Adipose",

    tolower(x) %in% c(
      "epithelial",
      "epithelium",
      "colon",
      "tumor",
      "tumour"
    ) ~ "Epithelial",

    TRUE ~ x
  )
}


read_cpdb_table <- function(path) {

  if (!file.exists(path)) {
    stop(
      "CellPhoneDB result file not found: ",
      path
    )
  }

  read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


find_pair_column <- function(
  data,
  sender,
  receiver
) {

  candidates <- c(
    paste0(
      sender,
      "|",
      receiver
    ),
    paste0(
      sender,
      ".",
      receiver
    ),
    paste0(
      sender,
      "_",
      receiver
    )
  )

  hit <- candidates[
    candidates %in% colnames(data)
  ]

  if (length(hit) == 0) {
    return(
      NA_character_
    )
  }

  hit[1]
}


###############################################################################
# 3. Load preprocessed GeoMx object
###############################################################################

if (!file.exists(PREPROCESSED_RDS)) {
  stop(
    "Preprocessed GeoMx object was not found:\n",
    PREPROCESSED_RDS,
    "\nRun 01_preprocessing_QC_normalization.R first."
  )
}

target_geomx <- readRDS(
  PREPROCESSED_RDS
)

if (!"q_norm" %in%
    Biobase::assayDataElementNames(
      target_geomx
    )) {
  stop(
    "The GeoMx object does not contain assay element 'q_norm'."
  )
}

metadata <- as.data.frame(
  Biobase::pData(
    target_geomx
  ),
  check.names = FALSE
)

message(
  "Loaded GeoMx object: ",
  nrow(target_geomx),
  " targets x ",
  ncol(target_geomx),
  " ROIs"
)


###############################################################################
# 4. Resolve metadata columns
###############################################################################

patient_col <- resolve_column(
  metadata,
  c(
    "Patient ID_from_CRC_base",
    "Patient.ID_from_CRC_base",
    "PatientID",
    "Patient_ID",
    "patient_id"
  ),
  "patient ID"
)

class_col <- resolve_column(
  metadata,
  c(
    "class",
    "Class"
  ),
  "normal/tumor class"
)

obesity_col <- resolve_column(
  metadata,
  c(
    "CategoryTwo",
    "Category",
    "Obesity",
    "obesity"
  ),
  "obesity group"
)

roi_type_col <- resolve_column(
  metadata,
  c(
    "ROIType",
    "ROI.Type",
    "roi_type",
    "segment"
  ),
  "ROI type"
)

segment_col <- if (
  "segment" %in%
    colnames(metadata)
) {
  "segment"
} else {
  NA_character_
}


###############################################################################
# 5. Build ROI-level CellPhoneDB metadata
###############################################################################

roi_info <- data.frame(
  barcode_original = clean_roi_id(
    rownames(metadata)
  ),
  PatientID = as.character(
    getElement(metadata, patient_col)
  ),
  Obesity = standardize_obesity(
    getElement(metadata, obesity_col)
  ),
  class = standardize_class(
    getElement(metadata, class_col)
  ),
  ROI_raw = as.character(
    getElement(metadata, roi_type_col)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!is.na(segment_col)) {

  roi_info$cell_type <- standardize_roi_type(
    roi_info$ROI_raw,
    getElement(metadata, segment_col)
  )

} else {

  roi_info$cell_type <- standardize_roi_type(
    roi_info$ROI_raw
  )
}

roi_info <- roi_info |>
  dplyr::filter(
    class == "Tumor",
    cell_type %in% c(
      "CD3",
      "Epithelial",
      "Adipose"
    ),
    !is.na(PatientID),
    PatientID != ""
  )


# Current CellPhoneDB documentation advises avoiding dashes in cell/barcode
# names. Keep the original ROI barcode in a mapping table, but use a safe
# internal barcode for CellPhoneDB.
roi_info <- roi_info |>
  dplyr::arrange(
    suppressWarnings(
      as.numeric(PatientID)
    ),
    cell_type,
    barcode_original
  ) |>
  dplyr::mutate(
    barcode_sample = paste0(
      "ROI_",
      sprintf(
        "%03d",
        dplyr::row_number()
      )
    )
  )


###############################################################################
# 6. Determine analyzable patients
###############################################################################

patient_compartments <- roi_info |>
  dplyr::group_by(
    PatientID
  ) |>
  dplyr::summarise(
    has_CD3 = any(
      cell_type == "CD3"
    ),
    has_Epithelial = any(
      cell_type == "Epithelial"
    ),
    has_Adipose = any(
      cell_type == "Adipose"
    ),
    n_ROIs = dplyr::n(),
    Obesity = dplyr::first(
      Obesity
    ),
    .groups = "drop"
  )

matched_patients <- patient_compartments |>
  dplyr::filter(
    has_CD3,
    has_Epithelial
  ) |>
  dplyr::pull(
    PatientID
  )


# Patients included in Figure 1H.

figure1h_patient_ids <- as.character(
  c(
    1, 3, 5, 6, 8, 9, 11, 13, 14, 15, 16, 18, 19, 20,
    22, 23, 25, 26, 28, 29, 30, 31, 32, 33, 34, 35, 36,
    39, 40, 41, 42, 44, 46, 49, 50, 51, 52, 53, 54, 56,
    57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85,
    86, 88, 89
  )
)

if (
  PATIENT_SET_MODE ==
    "figure1h"
) {

  selected_patients <- intersect(
    figure1h_patient_ids,
    matched_patients
  )

  missing_figure1h <- setdiff(
    figure1h_patient_ids,
    matched_patients
  )

  if (
    length(
      missing_figure1h
    ) > 0
  ) {
    warning(
      "The following Figure 1H patient IDs are not available as matched tumor CD3/Epithelial ROI pairs: ",
      paste(
        missing_figure1h,
        collapse = ", "
      )
    )
  }

} else if (
  PATIENT_SET_MODE ==
    "all_matched"
) {

  selected_patients <- matched_patients

} else {

  stop(
    "Unknown GEOMX_CPDB_PATIENT_SET_MODE: ",
    PATIENT_SET_MODE,
    ". Use 'figure1h' or 'all_matched'."
  )
}

if (
  length(
    selected_patients
  ) == 0
) {
  stop(
    "No matched tumor CD3/Epithelial patients were selected."
  )
}

roi_info_selected <- roi_info |>
  dplyr::filter(
    PatientID %in%
      selected_patients
  )

write.csv(
  roi_info,
  file.path(
    INPUT_DIR,
    "ROI_information_tumor_all.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  roi_info_selected,
  file.path(
    INPUT_DIR,
    "ROI_information_CellPhoneDB_selected.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  patient_compartments,
  file.path(
    INPUT_DIR,
    "patient_compartment_availability.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 7. Extract log2(Q3) expression matrix
###############################################################################

q3 <- Biobase::assayDataElement(
  target_geomx,
  elt = "q_norm"
)

q3 <- as.matrix(
  q3
)

storage.mode(q3) <- "numeric"

colnames(q3) <- clean_roi_id(
  colnames(q3)
)

feature_data <- as.data.frame(
  Biobase::fData(
    target_geomx
  ),
  check.names = FALSE
)

if (
  "CodeClass" %in%
    colnames(feature_data)
) {

  biological <- as.character(
    feature_data$CodeClass
  ) != "Negative"

} else {

  biological <- rep(
    TRUE,
    nrow(q3)
  )
}

gene_candidates <- c(
  "TargetName",
  "Gene",
  "gene",
  "GeneName",
  "gene_name"
)

gene_col <- gene_candidates[
  gene_candidates %in%
    colnames(feature_data)
][1]

if (
  !is.na(gene_col)
) {

  gene_names <- as.character(
    getElement(
      feature_data,
      gene_col
    )
  )

} else {

  gene_names <- rownames(q3)
}

q3_bio <- q3[
  biological,
  ,
  drop = FALSE
]

gene_names_bio <- gene_names[
  biological
]

if (
  anyNA(gene_names_bio) ||
    any(
      gene_names_bio == ""
    )
) {
  stop(
    "Missing biological gene names were detected."
  )
}

if (
  anyDuplicated(
    gene_names_bio
  )
) {
  stop(
    "Duplicated biological gene names were detected."
  )
}

rownames(q3_bio) <- gene_names_bio

if (
  any(
    q3_bio <= 0,
    na.rm = TRUE
  )
) {
  stop(
    "Non-positive Q3 values were detected."
  )
}

log2_q3 <- log2(
  q3_bio
)


###############################################################################
# 8. Match selected ROIs to expression columns
###############################################################################

missing_expression_rois <- setdiff(
  roi_info_selected$barcode_original,
  colnames(log2_q3)
)

if (
  length(
    missing_expression_rois
  ) > 0
) {
  stop(
    "Selected ROI(s) were not found in q_norm expression matrix: ",
    paste(
      missing_expression_rois,
      collapse = ", "
    )
  )
}


###############################################################################
# 9. Write CellPhoneDB inputs per patient
###############################################################################

for (
  patient_id in
    selected_patients
) {

  patient_meta <- roi_info_selected |>
    dplyr::filter(
      PatientID ==
        patient_id
    ) |>
    dplyr::arrange(
      factor(
        cell_type,
        levels = c(
          "CD3",
          "Epithelial",
          "Adipose"
        )
      ),
      barcode_original
    )

  ###########################################################################
  # CellPhoneDB metadata and expression input
  ###########################################################################

  patient_expr_original <- log2_q3[
    ,
    patient_meta$barcode_original,
    drop = FALSE
  ]

  ###########################################################################
  # CellPhoneDB-compatible format
  #
  # Use safe ROI barcodes in the CellPhoneDB input files.
  ###########################################################################

  cpdb_meta <- patient_meta |>
    dplyr::select(
      barcode_sample,
      cell_type
    )

  meta_file <- file.path(
    META_DIR,
    paste0(
      "ROI_",
      patient_id,
      ".txt"
    )
  )

  write.table(
    cpdb_meta,
    meta_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )


  patient_expr <- patient_expr_original

  colnames(
    patient_expr
  ) <- patient_meta$barcode_sample

  expr_out <- data.frame(
    Gene_Name = rownames(
      patient_expr
    ),
    as.data.frame(
      patient_expr,
      check.names = FALSE
    ),
    check.names = FALSE
  )

  expr_file <- file.path(
    EXPR_DIR,
    paste0(
      "expr_",
      patient_id,
      ".txt"
    )
  )

  write.table(
    expr_out,
    expr_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


write.csv(
  roi_info_selected |>
    dplyr::select(
      barcode_original,
      barcode_sample,
      PatientID,
      Obesity,
      class,
      cell_type
    ),
  file.path(
    INPUT_DIR,
    "ROI_barcode_mapping.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

message(
  "Prepared CellPhoneDB inputs for ",
  length(
    selected_patients
  ),
  " patients."
)


###############################################################################
# 10. Run CellPhoneDB Python worker
###############################################################################

if (
  RUN_CPDB
) {

  if (
    CPDB_FILE ==
      ""
  ) {
    stop(
      "GEOMX_CPDB_RUN=TRUE but GEOMX_CPDB_FILE was not set."
    )
  }

  if (
    !file.exists(
      CPDB_FILE
    )
  ) {
    stop(
      "CellPhoneDB database file was not found: ",
      CPDB_FILE
    )
  }

  if (
    !file.exists(
      PYTHON_SCRIPT
    )
  ) {
    stop(
      "CellPhoneDB Python worker was not found: ",
      PYTHON_SCRIPT
    )
  }

  python_args <- c(
    PYTHON_SCRIPT,
    "--cpdb-file",
    CPDB_FILE,
    "--meta-dir",
    META_DIR,
    "--expr-dir",
    EXPR_DIR,
    "--output-dir",
    RAW_DIR,
    "--iterations",
    as.character(
      CPDB_ITERATIONS
    ),
    "--threshold",
    as.character(
      CPDB_THRESHOLD
    ),
    "--threads",
    as.character(
      CPDB_THREADS
    ),
    "--counts-data",
    CPDB_COUNTS_DATA
  )

  message(
    "Running CellPhoneDB..."
  )

  exit_status <- system2(
    PYTHON_BIN,
    args = python_args
  )

  if (
    exit_status !=
      0
  ) {
    stop(
      "CellPhoneDB Python worker exited with status ",
      exit_status
    )
  }
}


###############################################################################
# 11. Read per-patient CellPhoneDB outputs, when available
###############################################################################

patient_output_dirs <- list.dirs(
  RAW_DIR,
  recursive = FALSE,
  full.names = TRUE
)

patient_output_dirs <- patient_output_dirs[
  grepl(
    "^Patient_",
    basename(
      patient_output_dirs
    )
  )
]

if (
  length(
    patient_output_dirs
  ) == 0
) {

  message("")
  message(
    "CellPhoneDB inputs were prepared successfully."
  )

  message(
    "No raw patient outputs were found yet under: ",
    RAW_DIR
  )

  message(
    "Run GeoMx/05_run_cellphonedb.py, or set GEOMX_CPDB_RUN=TRUE ",
    "and GEOMX_CPDB_FILE before rerunning this script."
  )

  stop(
    paste0(
      "CellPhoneDB inputs were prepared. ",
      "Run 05_run_cellphonedb.py (or set GEOMX_CPDB_RUN=TRUE), ",
      "then rerun 05_cellphonedb_analysis.R to merge outputs and draw Figure 1H."
    ),
    call. = FALSE
  )
}


###############################################################################
# 12. Build patient-direction long table
###############################################################################

roi_patient_group <- roi_info_selected |>
  dplyr::distinct(
    PatientID,
    Obesity
  )

direction_definitions <- tibble::tribble(
  ~Direction,        ~Sender,       ~Receiver,
  "Tumor - Immune",  "Epithelial",  "CD3",
  "Immune - Tumor",  "CD3",         "Epithelial"
)

patient_direction_results <- list()

result_counter <- 0L


for (
  patient_dir in
    patient_output_dirs
) {

  patient_id <- sub(
    "^Patient_",
    "",
    basename(
      patient_dir
    )
  )

  means_file <- file.path(
    patient_dir,
    "means.tsv"
  )

  pvalues_file <- file.path(
    patient_dir,
    "pvalues.tsv"
  )

  scores_file <- file.path(
    patient_dir,
    "interaction_scores.tsv"
  )

  if (
    !all(
      file.exists(
        c(
          means_file,
          pvalues_file,
          scores_file
        )
      )
    )
  ) {

    warning(
      "Skipping incomplete CellPhoneDB output for Patient ",
      patient_id
    )

    next
  }

  means_data <- read_cpdb_table(
    means_file
  )

  pvalues_data <- read_cpdb_table(
    pvalues_file
  )

  scores_data <- read_cpdb_table(
    scores_file
  )


  id_col <- if (
    "id_cp_interaction" %in%
      colnames(
        means_data
      )
  ) {
    "id_cp_interaction"
  } else if (
    "interaction_pair" %in%
      colnames(
        means_data
      )
  ) {
    "interaction_pair"
  } else {
    stop(
      "No CellPhoneDB interaction ID column found in ",
      means_file
    )
  }


  for (
    d in seq_len(
      nrow(
        direction_definitions
      )
    )
  ) {

    sender <- direction_definitions$Sender[
      d
    ]

    receiver <- direction_definitions$Receiver[
      d
    ]

    pair_col_mean <- find_pair_column(
      means_data,
      sender,
      receiver
    )

    pair_col_p <- find_pair_column(
      pvalues_data,
      sender,
      receiver
    )

    pair_col_score <- find_pair_column(
      scores_data,
      sender,
      receiver
    )

    if (
      anyNA(
        c(
          pair_col_mean,
          pair_col_p,
          pair_col_score
        )
      )
    ) {
      next
    }


    keep_metadata <- intersect(
      c(
        id_col,
        "interacting_pair",
        "partner_a",
        "partner_b",
        "gene_a",
        "gene_b",
        "secreted",
        "receptor_a",
        "receptor_b",
        "annotation_strategy",
        "is_integrin",
        "directionality",
        "classification"
      ),
      colnames(
        means_data
      )
    )


    current <- means_data |>
      dplyr::select(
        dplyr::all_of(
          keep_metadata
        )
      ) |>
      dplyr::mutate(
        MeanValue = as.numeric(
          getElement(means_data, pair_col_mean)
        ),
        PValue = as.numeric(
          getElement(pvalues_data, pair_col_p)
        ),
        Score = as.numeric(
          getElement(scores_data, pair_col_score)
        ),
        PatientID = patient_id,
        Direction = direction_definitions$Direction[
          d
        ],
        Sender = sender,
        Receiver = receiver
      )

    result_counter <- result_counter +
      1L

    patient_direction_results[result_counter] <- list(current)
  }
}


if (
  length(
    patient_direction_results
  ) == 0
) {
  stop(
    "No Epithelial/CD3 directional CellPhoneDB results were found."
  )
}


cpdb_long <- dplyr::bind_rows(
  patient_direction_results
) |>
  dplyr::left_join(
    roi_patient_group,
    by = "PatientID"
  ) |>
  dplyr::mutate(
    Obesity = factor(
      Obesity,
      levels = c(
        "Normal",
        "Obese"
      )
    )
  )


write.csv(
  cpdb_long,
  file.path(
    MERGED_DIR,
    "CellPhoneDB_patient_direction_long.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
###############################################################################

# A compact, robust wide export without relying on dplyr's dynamic rename.
wide_export <- function(
  data,
  direction,
  value_col,
  filename
) {

  sub <- data |>
    dplyr::filter(
      Direction ==
        direction
    )

  id_name <- if (
    "id_cp_interaction" %in%
      colnames(sub)
  ) {
    "id_cp_interaction"
  } else {
    "interaction_pair"
  }

  metadata_cols <- unique(
    c(
      id_name,
      intersect(
        c(
          "interacting_pair",
          "partner_a",
          "partner_b",
          "gene_a",
          "gene_b",
          "secreted",
          "receptor_a",
          "receptor_b",
          "annotation_strategy",
          "is_integrin",
          "directionality",
          "classification"
        ),
        colnames(sub)
      )
    )
  )

  meta_tbl <- sub |>
    dplyr::select(
      dplyr::all_of(
        metadata_cols
      )
    ) |>
    dplyr::distinct()

  value_tbl <- sub |>
    dplyr::transmute(
      interaction_id = !!rlang::sym(id_name),
      Patient = paste0(
        "Patient_",
        PatientID
      ),
      value = !!rlang::sym(value_col)
    ) |>
    tidyr::pivot_wider(
      names_from = Patient,
      values_from = value
    )

  names(
    meta_tbl
  )[
    names(
      meta_tbl
    ) ==
      id_name
  ] <- "interaction_id"

  out <- meta_tbl |>
    dplyr::left_join(
      value_tbl,
      by = "interaction_id"
    )

  names(out)[
    names(out) ==
      "interaction_id"
  ] <- "interaction_pair"

  write.table(
    out,
    file.path(
      MERGED_DIR,
      filename
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


for (
  direction_row in
    seq_len(
      nrow(
        direction_definitions
      )
    )
) {

  direction_name <- direction_definitions$Direction[
    direction_row
  ]

  prefix <- if (
    direction_name ==
      "Tumor - Immune"
  ) {
    "epithelial_CD3"
  } else {
    "CD3_epithelial"
  }

  wide_export(
    cpdb_long,
    direction_name,
    "MeanValue",
    paste0(
      prefix,
      "_means.txt"
    )
  )

  wide_export(
    cpdb_long,
    direction_name,
    "PValue",
    paste0(
      prefix,
      "_pvalues.txt"
    )
  )

  wide_export(
    cpdb_long,
    direction_name,
    "Score",
    paste0(
      prefix,
      "_interaction_scores.txt"
    )
  )
}


###############################################################################
# 14. Figure 1H data
###############################################################################

read_figure_table <- function(
  filename
) {

  path <- file.path(
    FIGURE1H_INPUT_DIR,
    filename
  )

  if (
    !file.exists(
      path
    )
  ) {
    stop(
      "Figure 1H input file not found: ",
      path
    )
  }

  x <- read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  good_names <- !is.na(
    colnames(x)
  ) &
    nzchar(
      colnames(x)
    )

  x[
    ,
    good_names,
    drop = FALSE
  ]
}


figure_patient_columns <- function(
  x
) {

  grep(
    "^Patient_[0-9]+$",
    colnames(x),
    value = TRUE
  )
}


prepare_figure_direction <- function(
  means_file,
  pvalues_file,
  scores_file,
  direction_label
) {

  means_data <- read_figure_table(
    means_file
  )

  pvalues_data <- read_figure_table(
    pvalues_file
  )

  scores_data <- read_figure_table(
    scores_file
  )

  patients <- figure_patient_columns(
    means_data
  )

  if (
    !identical(
      patients,
      figure_patient_columns(
        pvalues_data
      )
    ) ||
      !identical(
        patients,
        figure_patient_columns(
          scores_data
        )
      )
  ) {
    stop(
      "Patient columns differ among Figure 1H matrices."
    )
  }

  means_lr <- means_data |>
    dplyr::filter(
      directionality ==
        "Ligand-Receptor"
    ) |>
    dplyr::select(
      interacting_pair,
      dplyr::all_of(
        patients
      )
    )

  pvalues_lr <- pvalues_data |>
    dplyr::filter(
      directionality ==
        "Ligand-Receptor"
    ) |>
    dplyr::select(
      interacting_pair,
      dplyr::all_of(
        patients
      )
    )

  scores_lr <- scores_data |>
    dplyr::filter(
      directionality ==
        "Ligand-Receptor"
    ) |>
    dplyr::select(
      interacting_pair,
      dplyr::all_of(
        patients
      )
    )

  means_long <- means_lr |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patients
      ),
      names_to = "Patient",
      values_to = "MeanValue"
    )

  pvalues_long <- pvalues_lr |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patients
      ),
      names_to = "Patient",
      values_to = "PValue"
    )

  scores_long <- scores_lr |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        patients
      ),
      names_to = "Patient",
      values_to = "Score"
    )

  means_long |>
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
      Direction = direction_label
    )
}


figure_tumor_immune <- prepare_figure_direction(
  "epithelial_CD3_means.txt",
  "epithelial_CD3_pvalues.txt",
  "epithelial_CD3_interaction_scores.txt",
  "Tumor - Immune"
)

figure_immune_tumor <- prepare_figure_direction(
  "CD3_epithelial_means.txt",
  "CD3_epithelial_pvalues.txt",
  "CD3_epithelial_interaction_scores.txt",
  "Immune - Tumor"
)

figure_all <- dplyr::bind_rows(
  figure_tumor_immune,
  figure_immune_tumor
)


roi_annotation_path <- file.path(
  FIGURE1H_INPUT_DIR,
  "ROI_information_2.csv"
)

if (
  !file.exists(
    roi_annotation_path
  )
) {
  stop(
    "ROI_information_2.csv not found: ",
    roi_annotation_path
  )
}

figure_roi_info <- read.csv(
  roi_annotation_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

good_names <- !is.na(
  colnames(
    figure_roi_info
  )
) &
  nzchar(
    colnames(
      figure_roi_info
    )
  )

figure_roi_info <- figure_roi_info[
  ,
  good_names,
  drop = FALSE
]

figure_patient_info <- figure_roi_info |>
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

figure_all <- figure_all |>
  dplyr::left_join(
    figure_patient_info,
    by = "PatientID"
  )

figure1h_pairs <- tibble::tribble(
  ~Direction,       ~interacting_pair,              ~Display,
  "Tumor - Immune", "PVR_TIGIT",                    "PVR - TIGIT",
  "Tumor - Immune", "AREG_EGFR",                    "AREG - EGFR",
  "Tumor - Immune", "LTB_LTBR",                     "LTB - LTBR",
  "Tumor - Immune", "IGF1_IGF1R",                   "IGF1 - IGF1R",
  "Tumor - Immune", "EFNA4_EPHA1",                  "EFNA4 - EPHA1",
  "Immune - Tumor", "WNT8A_FZD5_LRP6",              "WNT8A - FZD5 - LRP6",
  "Immune - Tumor", "MICA_NKG2D_II_receptor",       "MICA - NKG2D II receptor",
  "Immune - Tumor", "JAG1_NOTCH3",                  "JAG1 - NOTCH3",
  "Immune - Tumor", "IL25_IL17_receptor_AB",        "IL25 - IL17 receptor AB",
  "Immune - Tumor", "GAS6_AXL",                     "GAS6 - AXL",
  "Immune - Tumor", "WNT10A_FZD5_LRP6",             "WNT10A - FZD5 - LRP6",
  "Immune - Tumor", "THBS2_CD36",                   "THBS2 - CD36",
  "Immune - Tumor", "IGF2_IGF1R",                   "IGF2 - IGF1R",
  "Immune - Tumor", "HLA-F_LILRB2",                 "HLA-F - LILRB2",
  "Immune - Tumor", "CTSG_FPR1",                    "CTSG - FPR1"
)

figure1h_data <- figure_all |>
  dplyr::inner_join(
    figure1h_pairs,
    by = c(
      "Direction",
      "interacting_pair"
    )
  ) |>
  dplyr::filter(
    PValue ==
      0
  )

figure_patient_order <- figure_all |>
  dplyr::distinct(
    PatientID,
    Obesity
  ) |>
  dplyr::mutate(
    PatientNumeric = suppressWarnings(
      as.numeric(
        PatientID
      )
    ),
    GroupOrder = ifelse(
      Obesity ==
        "Non-Obese",
      1,
      2
    )
  ) |>
  dplyr::arrange(
    GroupOrder,
    PatientNumeric
  ) |>
  dplyr::pull(
    PatientID
  )

figure1h_data <- figure1h_data |>
  dplyr::mutate(
    PatientID = factor(
      PatientID,
      levels =
        figure_patient_order
    ),
    Obesity = factor(
      Obesity,
      levels = c(
        "Non-Obese",
        "Obese"
      )
    )
  )

write.csv(
  figure1h_data,
  file.path(
    MERGED_DIR,
    "Figure1H_plot_data.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


###############################################################################
# 15. Figure 1H
###############################################################################

tumor_immune_levels <- figure1h_pairs |>
  dplyr::filter(
    Direction ==
      "Tumor - Immune"
  ) |>
  dplyr::pull(
    Display
  )

immune_tumor_levels <- figure1h_pairs |>
  dplyr::filter(
    Direction ==
      "Immune - Tumor"
  ) |>
  dplyr::pull(
    Display
  )


make_direction_plot <- function(
  data,
  direction,
  interaction_levels
) {

  plot_data <- data |>
    dplyr::filter(
      Direction ==
        direction
    ) |>
    dplyr::mutate(
      Display = factor(
        Display,
        levels = rev(
          interaction_levels
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
      ),
      alpha = 0.95
    ) +
    facet_grid(
      . ~ Obesity,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_color_gradient(
      low = "lightblue",
      high = "red",
      name = "Mean"
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
        4.8
      ),
      name = "Interaction\nScore"
    ) +
    labs(
      title = direction,
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      axis.text.x =
        element_blank(),
      axis.ticks.x =
        element_blank(),
      axis.text.y =
        element_text(
          size = 10
        ),
      strip.background =
        element_rect(
          fill = "white",
          color = "black",
          linewidth = 0.5
        ),
      strip.text =
        element_text(
          size = 10
        ),
      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold",
          size = 13
        ),
      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.6
        ),
      panel.grid.major =
        element_line(
          color = "grey92",
          linewidth = 0.25
        ),
      panel.grid.minor =
        element_blank(),
      legend.position =
        "right"
    )
}


plot_tumor_immune <- make_direction_plot(
  figure1h_data,
  "Tumor - Immune",
  tumor_immune_levels
)

plot_immune_tumor <- make_direction_plot(
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
  file.path(
    FIGURE_DIR,
    "Figure1H_CellPhoneDB.pdf"
  ),
  figure1h,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  file.path(
    FIGURE_DIR,
    "Figure1H_CellPhoneDB.png"
  ),
  figure1h,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)


###############################################################################
# 16. Summaries
###############################################################################

patient_summary <- roi_info_selected |>
  dplyr::distinct(
    PatientID,
    Obesity
  ) |>
  dplyr::count(
    Obesity,
    name = "N_patients"
  )

write.csv(
  patient_summary,
  file.path(
    OUTPUT_DIR,
    "CellPhoneDB_patient_counts.csv"
  ),
  row.names = FALSE
)

summary_lines <- c(
  paste0(
    "Input GeoMx object: ",
    PREPROCESSED_RDS
  ),
  paste0(
    "Patient-set mode: ",
    PATIENT_SET_MODE
  ),
  paste0(
    "Patients selected: ",
    length(
      selected_patients
    )
  ),
  paste0(
    "Tumor ROIs supplied to CellPhoneDB: ",
    nrow(
      roi_info_selected
    )
  ),
  "Representation: one GeoMx ROI = one CellPhoneDB observation/pseudo-cell",
  "CellPhoneDB cluster labels: CD3, Epithelial, Adipose",
  "No pseudobulk averaging across ROIs was performed",
  "Expression input: log2(Q3-normalized expression)",
  paste0(
    "CellPhoneDB iterations: ",
    CPDB_ITERATIONS
  ),
  paste0(
    "CellPhoneDB threshold: ",
    CPDB_THRESHOLD
  ),
  paste0(
    "CellPhoneDB counts_data: ",
    CPDB_COUNTS_DATA
  ),
  "Figure 1H significance filter: PValue == 0",
  "Tumor - Immune direction: Epithelial -> CD3",
  "Immune - Tumor direction: CD3 -> Epithelial"
)

writeLines(
  summary_lines,
  file.path(
    OUTPUT_DIR,
    "CellPhoneDB_analysis_summary.txt"
  )
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


message("")
message(
  "GeoMx CellPhoneDB analysis / Figure 1H completed."
)

message(
  "Results written to: ",
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  )
)
