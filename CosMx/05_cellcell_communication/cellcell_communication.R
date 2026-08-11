#!/usr/bin/env Rscript

###############################################################################
# 05_CosMx_cellcell_communication.R
#
# CosMx cell-cell communication analysis using CellChat
# Obesity-associated colorectal cancer study
#
#
# Workflow
# --------
# 1. Load annotated CosMx Seurat object
# 2. Set final cell-type identities
# 3. Split patients into obese and non-obese groups
# 4. Run spatial CellChat separately for each patient
# 5. Compute pathway-level communication and aggregate networks
# 6. Calculate group-average communication networks
# 7. Compare obese vs non-obese interaction strength
# 8. Compare Tumor -> SPP1+ pathway-level communication
# 9. Save CellChat objects and summary tables
#
# Individual ligand-receptor pair analysis is performed separately in:
#   06_CosMx_ligand_receptor_analysis.R
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- "path/to/CosMx_annotated_object.rds"

OUTPUT_DIR <- "results/CosMx/05_cellcell_communication"

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
library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(future)


###############################################################################
# 2. Load annotated CosMx Seurat object
###############################################################################

cosmx <- readRDS(
  INPUT_RDS
)


###############################################################################
# 3. Metadata and analysis settings
###############################################################################

# Detailed cell-type annotation used in the original CellChat workflow.
CELLTYPE_COLUMN <- "cell_type_figure5_032825"

PATIENT_COLUMN <- "patient"

OBESITY_COLUMN <- "obese_status"


# Spatial CellChat parameters from the original analysis.
CONVERSION_FACTOR <- 0.18

INTERACTION_RANGE <- 250

SCALE_DISTANCE <- 575

CONTACT_RANGE <- 20


# Minimum number of valid cells required to run one patient.
MIN_CELLS <- 100


# Original CellChat probability settings.
COMMUNICATION_TYPE <- "truncatedMean"

TRIM <- 0.1


###############################################################################
# 4. Set cell-type identities
#
# The original CellChat script used Idents(cosmx) as the CellChat label.
# This assignment is made explicit here so that the script does not depend
# on a previously stored active identity.
###############################################################################

Idents(cosmx) <- CELLTYPE_COLUMN


###############################################################################
# 5. CellChat database and execution mode
###############################################################################

CellChatDB.use <- subsetDB(
  CellChatDB.human
)


# Original analysis was run sequentially.
future::plan(
  "sequential"
)


###############################################################################
# 6. Identify obese and non-obese patients
###############################################################################

meta <- cosmx@meta.data


patient_status <- meta %>%
  filter(
    !is.na(.data[[PATIENT_COLUMN]]),
    !is.na(.data[[OBESITY_COLUMN]])
  ) %>%
  distinct(
    .data[[PATIENT_COLUMN]],
    .data[[OBESITY_COLUMN]]
  )


obese_patients <- patient_status %>%
  filter(
    .data[[OBESITY_COLUMN]] == "obese"
  ) %>%
  pull(
    all_of(PATIENT_COLUMN)
  )


nonobese_patients <- patient_status %>%
  filter(
    .data[[OBESITY_COLUMN]] == "non-obese"
  ) %>%
  pull(
    all_of(PATIENT_COLUMN)
  )


###############################################################################
# 7. Function: run spatial CellChat for one patient
#
# Original workflow:
#
#   data.input = SCT normalized expression
#   coordinates = x_slide_mm / y_slide_mm
#   datatype = "spatial"
#
#   identifyOverExpressedGenes()
#   identifyOverExpressedInteractions(variable.both = FALSE)
#
#   computeCommunProb(
#       type = "truncatedMean",
#       trim = 0.1,
#       distance.use = TRUE,
#       interaction.range = 250,
#       scale.distance = 575,
#       contact.dependent = TRUE,
#       contact.range = 20
#   )
###############################################################################

run_cellchat_patient <- function(
    object,
    patient_id
) {

  cells_patient <- rownames(
    object@meta.data
  )[
    object@meta.data[[PATIENT_COLUMN]] == patient_id
  ]


  valid_cells <- intersect(
    cells_patient,
    colnames(
      GetAssayData(
        object,
        slot = "data",
        assay = "SCT"
      )
    )
  )


  if (
    length(valid_cells) < MIN_CELLS
  ) {

    message(
      "Skipping ",
      patient_id,
      " (too few valid cells)"
    )

    return(NULL)
  }


  data.input <- GetAssayData(
    object,
    slot = "data",
    assay = "SCT"
  )[
    ,
    valid_cells,
    drop = FALSE
  ]


  meta.input <- data.frame(
    labels = Idents(object)[valid_cells],
    samples = patient_id,
    row.names = valid_cells
  )


  coords <- object@meta.data[
    valid_cells,
    c(
      "x_slide_mm",
      "y_slide_mm"
    )
  ]


  colnames(coords) <- c(
    "x",
    "y"
  )


  d <- computeCellDistance(
    coords
  )


  spot.size <- min(
    d[d > 0],
    na.rm = TRUE
  ) * CONVERSION_FACTOR


  spatial.factors <- data.frame(
    ratio = CONVERSION_FACTOR,
    tol = spot.size / 2
  )


  cellchat <- createCellChat(
    object = data.input,
    meta = meta.input,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = spatial.factors
  )


  cellchat@DB <- CellChatDB.use


  cellchat <- subsetData(
    cellchat
  )


  cellchat <- identifyOverExpressedGenes(
    cellchat
  )


  cellchat <- identifyOverExpressedInteractions(
    cellchat,
    variable.both = FALSE
  )


  cellchat <- computeCommunProb(
    cellchat,
    type = COMMUNICATION_TYPE,
    trim = TRIM,
    distance.use = TRUE,
    interaction.range = INTERACTION_RANGE,
    scale.distance = SCALE_DISTANCE,
    contact.dependent = TRUE,
    contact.range = CONTACT_RANGE
  )


  cellchat <- computeCommunProbPathway(
    cellchat
  )


  cellchat <- aggregateNet(
    cellchat
  )


  return(
    cellchat
  )
}


###############################################################################
# 8. Run CellChat separately for each obese patient
###############################################################################

cellchat_list_obese <- list()


for (
  pt in obese_patients
) {

  message(
    "Processing obese patient: ",
    pt
  )


  result <- run_cellchat_patient(
    object = cosmx,
    patient_id = pt
  )


  if (
    !is.null(result)
  ) {

    cellchat_list_obese[[pt]] <- result
  }
}


###############################################################################
# 9. Run CellChat separately for each non-obese patient
###############################################################################

cellchat_list_nonobese <- list()


for (
  pt in nonobese_patients
) {

  message(
    "Processing non-obese patient: ",
    pt
  )


  result <- run_cellchat_patient(
    object = cosmx,
    patient_id = pt
  )


  if (
    !is.null(result)
  ) {

    cellchat_list_nonobese[[pt]] <- result
  }
}


###############################################################################
# 10. Save patient-level CellChat objects
###############################################################################

saveRDS(
  cellchat_list_obese,
  file.path(
    OUTPUT_DIR,
    "CellChat_obese_patient_list.rds"
  )
)


saveRDS(
  cellchat_list_nonobese,
  file.path(
    OUTPUT_DIR,
    "CellChat_nonobese_patient_list.rds"
  )
)


###############################################################################
# 11. Function: calculate mean interaction-weight network
#
# The original workflow extracted @net$weight from each patient and averaged
# the matrices across patients after restricting to common cell types.
###############################################################################

get_mean_weight_network <- function(
    cellchat_list
) {

  net_list <- lapply(
    cellchat_list,
    function(x) {
      x@net$weight
    }
  )


  net_list <- Filter(
    function(x) {
      !is.null(x) &&
        nrow(x) > 0
    },
    net_list
  )


  if (
    length(net_list) == 0
  ) {

    stop(
      "No valid CellChat interaction-weight matrices were found."
    )
  }


  common_celltypes <- Reduce(
    intersect,
    lapply(
      net_list,
      rownames
    )
  )


  net_list <- lapply(
    net_list,
    function(x) {

      x[
        common_celltypes,
        common_celltypes,
        drop = FALSE
      ]
    }
  )


  mean_net <- Reduce(
    "+",
    net_list
  ) / length(
    net_list
  )


  return(
    mean_net
  )
}


###############################################################################
# 12. Group-average communication networks
###############################################################################

mean_net_obese <- get_mean_weight_network(
  cellchat_list_obese
)


mean_net_nonobese <- get_mean_weight_network(
  cellchat_list_nonobese
)


###############################################################################
# 13. Align cell types between obese and non-obese networks
###############################################################################

common_celltypes <- intersect(
  rownames(mean_net_obese),
  rownames(mean_net_nonobese)
)


mean_net_obese <- mean_net_obese[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]


mean_net_nonobese <- mean_net_nonobese[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]


###############################################################################
# 14. Save group-average interaction matrices
###############################################################################

write.csv(
  mean_net_obese,
  file.path(
    OUTPUT_DIR,
    "mean_interaction_weight_obese.csv"
  )
)


write.csv(
  mean_net_nonobese,
  file.path(
    OUTPUT_DIR,
    "mean_interaction_weight_nonobese.csv"
  )
)


###############################################################################
# 15. Obese and non-obese global communication networks
###############################################################################

shared_max <- max(
  mean_net_obese,
  mean_net_nonobese,
  na.rm = TRUE
)


jpeg(
  file.path(
    OUTPUT_DIR,
    "plots",
    "CellChat_network_obese.jpeg"
  ),
  width = 1800,
  height = 1800,
  res = 300
)


netVisual_circle(
  mean_net_obese,
  vertex.weight = rowSums(
    mean_net_obese
  ),
  weight.scale = TRUE,
  edge.weight.max = shared_max,
  label.edge = FALSE,
  title.name = "Average interaction weight: Obese"
)


dev.off()


jpeg(
  file.path(
    OUTPUT_DIR,
    "plots",
    "CellChat_network_nonobese.jpeg"
  ),
  width = 1800,
  height = 1800,
  res = 300
)


netVisual_circle(
  mean_net_nonobese,
  vertex.weight = rowSums(
    mean_net_nonobese
  ),
  weight.scale = TRUE,
  edge.weight.max = shared_max,
  label.edge = FALSE,
  title.name = "Average interaction weight: Non-obese"
)


dev.off()


###############################################################################
# 16. Network difference: Obese - Non-obese
###############################################################################

diff_weight <- mean_net_obese -
  mean_net_nonobese


diff_weight[
  !is.finite(
    diff_weight
  )
] <- 0


diff_weight[
  abs(
    diff_weight
  ) < 1e-5
] <- 0


write.csv(
  diff_weight,
  file.path(
    OUTPUT_DIR,
    "interaction_weight_difference_obese_minus_nonobese.csv"
  )
)


jpeg(
  file.path(
    OUTPUT_DIR,
    "plots",
    "CellChat_network_difference_heatmap.jpeg"
  ),
  width = 1800,
  height = 1600,
  res = 300
)


pheatmap::pheatmap(
  diff_weight,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Interaction Weight Difference (Obese - Non-obese)"
)


dev.off()


###############################################################################
# 17. Tumor <-> SPP1+ network
###############################################################################

key_cells <- c(
  "Tumor Cells",
  "SPP1+"
)


key_cells <- intersect(
  key_cells,
  common_celltypes
)


if (
  length(key_cells) == 2
) {

  tumor_spp1_obese <- mean_net_obese[
    key_cells,
    key_cells,
    drop = FALSE
  ]


  tumor_spp1_nonobese <- mean_net_nonobese[
    key_cells,
    key_cells,
    drop = FALSE
  ]


  shared_max_tumor_spp1 <- max(
    tumor_spp1_obese,
    tumor_spp1_nonobese,
    na.rm = TRUE
  )


  jpeg(
    file.path(
      OUTPUT_DIR,
      "plots",
      "Tumor_SPP1_network_obese.jpeg"
    ),
    width = 1400,
    height = 1400,
    res = 300
  )


  netVisual_circle(
    tumor_spp1_obese,
    vertex.weight = rowSums(
      tumor_spp1_obese
    ),
    weight.scale = TRUE,
    edge.weight.max = shared_max_tumor_spp1,
    label.edge = FALSE,
    title.name = "Obese: Tumor - SPP1+"
  )


  dev.off()


  jpeg(
    file.path(
      OUTPUT_DIR,
      "plots",
      "Tumor_SPP1_network_nonobese.jpeg"
    ),
    width = 1400,
    height = 1400,
    res = 300
  )


  netVisual_circle(
    tumor_spp1_nonobese,
    vertex.weight = rowSums(
      tumor_spp1_nonobese
    ),
    weight.scale = TRUE,
    edge.weight.max = shared_max_tumor_spp1,
    label.edge = FALSE,
    title.name = "Non-obese: Tumor - SPP1+"
  )


  dev.off()
}


###############################################################################
# 18. Function: average pathway-level communication across patients
#
# Original workflow:
#
#   subsetCommunication(object, slot.name = "netP")
#   aggregate(prob ~ pathway_name + source + target, FUN = mean)
###############################################################################

get_mean_pathway_network <- function(
    cellchat_list
) {

  pathway_tables <- lapply(
    cellchat_list,
    function(x) {

      subsetCommunication(
        x,
        slot.name = "netP"
      )
    }
  )


  pathway_tables <- Filter(
    function(x) {
      !is.null(x) &&
        nrow(x) > 0
    },
    pathway_tables
  )


  pathway_data <- bind_rows(
    pathway_tables
  )


  pathway_mean <- pathway_data %>%
    group_by(
      pathway_name,
      source,
      target
    ) %>%
    summarise(
      prob = mean(
        prob,
        na.rm = TRUE
      ),
      .groups = "drop"
    )


  return(
    pathway_mean
  )
}


###############################################################################
# 19. Obese and non-obese pathway-level networks
###############################################################################

pathway_obese <- get_mean_pathway_network(
  cellchat_list_obese
)


pathway_nonobese <- get_mean_pathway_network(
  cellchat_list_nonobese
)


write.csv(
  pathway_obese,
  file.path(
    OUTPUT_DIR,
    "pathway_communication_obese.csv"
  ),
  row.names = FALSE
)


write.csv(
  pathway_nonobese,
  file.path(
    OUTPUT_DIR,
    "pathway_communication_nonobese.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 20. Tumor -> SPP1+ pathway-level comparison
###############################################################################

tumor_spp1_pathway_obese <- pathway_obese %>%
  filter(
    source == "Tumor Cells",
    target == "SPP1+"
  )


tumor_spp1_pathway_nonobese <- pathway_nonobese %>%
  filter(
    source == "Tumor Cells",
    target == "SPP1+"
  )


pathway_comparison <- merge(
  tumor_spp1_pathway_obese,
  tumor_spp1_pathway_nonobese,
  by = "pathway_name",
  suffixes = c(
    "_obese",
    "_nonobese"
  )
)


pathway_comparison <- pathway_comparison %>%
  mutate(
    log2FC = log2(
      (
        prob_obese +
          1e-6
      ) /
        (
          prob_nonobese +
            1e-6
        )
    )
  ) %>%
  arrange(
    desc(
      log2FC
    )
  )


write.csv(
  pathway_comparison,
  file.path(
    OUTPUT_DIR,
    "Tumor_to_SPP1_pathway_comparison.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 21. Tumor -> SPP1+ pathway comparison plot
#
# Above the diagonal:
#   stronger pathway-level communication in obese tumors.
#
# Below the diagonal:
#   stronger pathway-level communication in non-obese tumors.
###############################################################################

if (
  nrow(pathway_comparison) > 0
) {

  p_pathway <- ggplot(
    pathway_comparison,
    aes(
      x = prob_nonobese,
      y = prob_obese
    )
  ) +
    geom_point(
      aes(
        size = abs(
          log2FC
        )
      ),
      alpha = 0.8
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    geom_text_repel(
      aes(
        label = pathway_name
      ),
      size = 3.5,
      max.overlaps = 30
    ) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Communication Probability (Non-obese)",
      y = "Communication Probability (Obese)",
      size = "|log2FC|",
      title = "Tumor -> SPP1+ Pathway Communication"
    ) +
    theme_minimal()


  ggsave(
    file.path(
      OUTPUT_DIR,
      "plots",
      "Tumor_to_SPP1_pathway_obese_vs_nonobese.jpeg"
    ),
    plot = p_pathway,
    width = 7,
    height = 6,
    dpi = 300
  )
}


###############################################################################
# 22. Save combined CellChat results
###############################################################################

cellchat_results <- list(
  obese = cellchat_list_obese,
  nonobese = cellchat_list_nonobese,
  mean_net_obese = mean_net_obese,
  mean_net_nonobese = mean_net_nonobese,
  diff_weight = diff_weight,
  pathway_obese = pathway_obese,
  pathway_nonobese = pathway_nonobese,
  tumor_to_SPP1_pathway = pathway_comparison
)


saveRDS(
  cellchat_results,
  file.path(
    OUTPUT_DIR,
    "CosMx_CellChat_obese_vs_nonobese.rds"
  )
)
