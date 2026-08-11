#!/usr/bin/env Rscript

###############################################################################
# 06_CosMx_ligand_receptor_analysis.R
#
# CosMx ligand-receptor analysis
#
# IMPORTANT
# ---------
# SPP1+ is a refined macrophage subtype.
# Therefore this script uses only CellChat objects generated with:
#
#   cell_type_refined
#
# from:
#   05_CosMx_cellcell_communication.R
#
# Workflow
# --------
# 1. Load refined patient-level CellChat results
# 2. Extract ligand-receptor interactions for each patient
# 3. Select Tumor Cells -> SPP1+ interactions
# 4. Compare communication probability in Obese vs Non-obese
# 5. Calculate interaction-level log2 fold change
# 6. Perform patient-level Wilcoxon tests
# 7. Save differential ligand-receptor tables
###############################################################################


###############################################################################
# 0. User settings
###############################################################################

INPUT_RDS <- paste0(
  "results/CosMx/05_cellcell_communication/",
  "CosMx_CellChat_obese_vs_nonobese.rds"
)

OUTPUT_DIR <- "results/CosMx/06_ligand_receptor_analysis"

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

library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)


###############################################################################
# 2. Load CellChat results
###############################################################################

cellchat_results <- readRDS(
  INPUT_RDS
)

cellchat_obese <- cellchat_results$refined$obese

cellchat_nonobese <- cellchat_results$refined$nonobese


###############################################################################
# 3. Extract LR communication table from each patient
###############################################################################

extract_lr_patient <- function(
    cellchat_list,
    obesity_group
) {

  result <- lapply(
    names(cellchat_list),
    function(pt) {

      df <- subsetCommunication(
        cellchat_list[[pt]],
        slot.name = "net"
      )

      if (
        is.null(df) ||
        nrow(df) == 0
      ) {
        return(NULL)
      }

      df$patient <- pt
      df$obese_status <- obesity_group

      df
    }
  )

  result <- Filter(
    Negate(is.null),
    result
  )

  bind_rows(
    result
  )
}


lr_obese <- extract_lr_patient(
  cellchat_obese,
  "obese"
)

lr_nonobese <- extract_lr_patient(
  cellchat_nonobese,
  "non-obese"
)

lr_all <- bind_rows(
  lr_obese,
  lr_nonobese
)


write.csv(
  lr_all,
  file.path(
    OUTPUT_DIR,
    "all_refined_ligand_receptor_interactions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 4. Tumor Cells -> SPP1+ ligand-receptor interactions
###############################################################################

tumor_spp1_lr <- lr_all %>%
  filter(
    source == "Tumor Cells",
    target == "SPP1+"
  )


write.csv(
  tumor_spp1_lr,
  file.path(
    OUTPUT_DIR,
    "Tumor_to_SPP1_patient_level_interactions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 5. Interaction identifier
#
# CellChat normally provides interaction_name and interaction_name_2.
# interaction_name_2 is preferred for readable plotting when available.
###############################################################################

if (
  "interaction_name_2" %in%
    colnames(tumor_spp1_lr)
) {

  tumor_spp1_lr$interaction_id <-
    tumor_spp1_lr$interaction_name_2

} else {

  tumor_spp1_lr$interaction_id <-
    tumor_spp1_lr$interaction_name
}


###############################################################################
# 6. Mean communication probability by obesity group
###############################################################################

lr_group_mean <- tumor_spp1_lr %>%
  group_by(
    interaction_id,
    obese_status
  ) %>%
  summarise(
    mean_prob = mean(
      prob,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = obese_status,
    values_from = mean_prob,
    values_fill = 0
  ) %>%
  rename(
    prob_obese = obese,
    prob_nonobese = `non-obese`
  ) %>%
  mutate(
    log2FC = log2(
      (prob_obese + 1e-6) /
        (prob_nonobese + 1e-6)
    )
  )


###############################################################################
# 7. Patient-level Wilcoxon test for each LR interaction
#
# Missing interaction in a patient is treated as communication probability 0.
###############################################################################

patients_obese <- names(
  cellchat_obese
)

patients_nonobese <- names(
  cellchat_nonobese
)

all_patients <- c(
  patients_obese,
  patients_nonobese
)

patient_info <- data.frame(
  patient = all_patients,
  obese_status = c(
    rep(
      "obese",
      length(patients_obese)
    ),
    rep(
      "non-obese",
      length(patients_nonobese)
    )
  ),
  stringsAsFactors = FALSE
)


lr_patient_prob <- tumor_spp1_lr %>%
  group_by(
    patient,
    obese_status,
    interaction_id
  ) %>%
  summarise(
    prob = mean(
      prob,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


all_interactions <- unique(
  lr_patient_prob$interaction_id
)


lr_complete <- expand_grid(
  patient = patient_info$patient,
  interaction_id = all_interactions
) %>%
  left_join(
    patient_info,
    by = "patient"
  ) %>%
  left_join(
    lr_patient_prob,
    by = c(
      "patient",
      "obese_status",
      "interaction_id"
    )
  ) %>%
  mutate(
    prob = replace_na(
      prob,
      0
    )
  )


lr_pvalues <- lr_complete %>%
  group_by(
    interaction_id
  ) %>%
  group_modify(
    ~ {

      dat <- .x

      if (
        length(
          unique(
            dat$obese_status
          )
        ) < 2
      ) {

        return(
          data.frame(
            p_value = NA_real_
          )
        )
      }

      data.frame(
        p_value = wilcox.test(
          prob ~ obese_status,
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
    ),
    minus_log10_p = -log10(
      pmax(
        p_value,
        1e-300
      )
    )
  )


###############################################################################
# 8. Final differential LR table
###############################################################################

lr_differential <- lr_group_mean %>%
  left_join(
    lr_pvalues,
    by = "interaction_id"
  ) %>%
  arrange(
    desc(
      log2FC
    )
  )


write.csv(
  lr_differential,
  file.path(
    OUTPUT_DIR,
    "Tumor_to_SPP1_LR_obese_vs_nonobese.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 9. LR log2FC plot
#
# Positive log2FC:
#   stronger in Obese
#
# Negative log2FC:
#   stronger in Non-obese
###############################################################################

plot_lr <- lr_differential %>%
  filter(
    is.finite(log2FC)
  ) %>%
  arrange(
    log2FC
  ) %>%
  mutate(
    interaction_id = factor(
      interaction_id,
      levels = interaction_id
    )
  )


p_fc <- ggplot(
  plot_lr,
  aes(
    x = interaction_id,
    y = log2FC
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "",
    y = "log2FC (Obese / Non-obese)",
    title = "Tumor Cells -> SPP1+ ligand-receptor interactions"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "Tumor_to_SPP1_LR_log2FC.jpeg"
  ),
  p_fc,
  width = 12,
  height = 5,
  dpi = 300
)


###############################################################################
# 10. Communication-probability heatmap-style plot
###############################################################################

lr_heatmap_data <- lr_group_mean %>%
  select(
    interaction_id,
    prob_obese,
    prob_nonobese
  ) %>%
  pivot_longer(
    cols = c(
      prob_obese,
      prob_nonobese
    ),
    names_to = "group",
    values_to = "prob"
  ) %>%
  mutate(
    group = recode(
      group,
      prob_nonobese = "Non-obese",
      prob_obese = "Obese"
    )
  )


p_heat <- ggplot(
  lr_heatmap_data,
  aes(
    x = group,
    y = interaction_id,
    fill = prob
  )
) +
  geom_tile() +
  labs(
    x = "",
    y = "",
    fill = "Comm. Prob"
  ) +
  theme_minimal()


ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    "Tumor_to_SPP1_LR_communication_probability.jpeg"
  ),
  p_heat,
  width = 5,
  height = 10,
  dpi = 300
)


###############################################################################
# 11. Save analysis object
###############################################################################

lr_results <- list(
  refined_annotation = "cell_type_refined",
  patient_level = tumor_spp1_lr,
  patient_complete = lr_complete,
  differential = lr_differential
)

saveRDS(
  lr_results,
  file.path(
    OUTPUT_DIR,
    "CosMx_Tumor_to_SPP1_LR_results.rds"
  )
)
