# ==============================================================================
# 16_summarise_glioma_sensitivity.R
# Run after files 14 and 15: Rscript 16_summarise_glioma_sensitivity.R
# Publication tables and figures for target-anchored glioma
# VITA-HDM sensitivity analyses.
# Inputs:
#   results/case_study/sensitivity/tables/
# Outputs:
#   results/case_study/sensitivity/final_outputs/
# Output specification:
#   - Keep all final sensitivity tables for SA1-SA8.
#   - Generate a concise 4-panel main-text sensitivity figure:
#       A. Noise knock-in robustness table
#       B. Source dose-response gene-set recovery
#       C. BFDR/stability threshold grid
#       D. Repeated-CV predictive gain of fixed signatures
#   - Generate separate high-resolution files for each main-text panel.
#   - Generate supplementary figures for detailed diagnostics:
#       SA2 C-index dose-response
#       SA3 BFDR/stability threshold grid
#       SA5 Cox-Snell residual diagnostics
#       SA6 raw source-target PCA
#       SA6 within-cohort standardized PCA
#       SA7 hyperparameter grid
#       optional 7-panel sensitivity overview
# Dependencies:
#   - Required R packages: ggplot2, dplyr, tidyr, scales, patchwork.
#   - Optional R package: ggrepel (plain-text fallback is used if unavailable).
#   - Required project scripts: 00_project_setup.R,
#     08_glioma_vita_functions.R, and 13_glioma_sensitivity_common.R.
# Notes:
#   - This script does not refit models and does not rerun bootstrap analyses.
#   - It reads existing SA1-SA8 tables, refreshes baseline classification
#     tables through the shared context, and regenerates final outputs.
#   - Explicit namespaces are used to avoid select()/filter()/mutate() conflicts.
#   - No sink() logging is used.
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
# ------------------------------------------------------------------------------
# Project paths
# ------------------------------------------------------------------------------
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(
    sub("^--file=", "", script_file_arg[[1]]),
    winslash = "/",
    mustWork = FALSE
  ))
} else if (!is.null(source_file)) {
  dirname(normalizePath(
    source_file,
    winslash = "/",
    mustWork = FALSE
  ))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(SCRIPT_DIR, "00_project_setup.R"))
PROJECT_ROOT <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_ROOT)
require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist:\n  ", path, call. = FALSE)
  }
  normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
}
FINAL_SELECTION_RDS_ABS <- require_file(
  file.path(
    PATHS$discovery,
    "glioma_final_selection_objects.rds"
  ),
  "Final selection object"
)
TARGET_ANCHORED_REFERENCE_FILE <- require_file(
  file.path(
    PATHS$discovery,
    "glioma_final_selected_genes_target_anchored.csv"
  ),
  "Target-anchored reference table"
)
COMMON_SCRIPT_PATH <- require_file(
  file.path(
    PROJECT_ROOT,
    "13_glioma_sensitivity_common.R"
  ),
  "Common sensitivity script"
)
require_file(
  file.path(
    PROJECT_ROOT,
    "08_glioma_vita_functions.R"
  ),
  "VITA-HDM function script"
)
if (!dir.exists(PATHS$sensitivity_tables)) {
  stop(
    "Sensitivity input table directory does not exist:\n  ",
    PATHS$sensitivity_tables,
    "\nRun files 14 and 15 first.",
    call. = FALSE
  )
}
source(
  COMMON_SCRIPT_PATH,
  local = .GlobalEnv,
  chdir = FALSE
)
RESULTS_DIR <- PATHS$case_results
SENS_DIR <- PATHS$sensitivity
SENS_TAB_DIR <- PATHS$sensitivity_tables
SENS_PRIVATE_DIR <- PATHS$sensitivity_private
SENS_FINAL_DIR <- ensure_directory(PATHS$sensitivity_final)
# Load the canonical seven-gene order directly from the absolute discovery file.
# This order is used only for display; gene labels are never converted to a
# factor, so a mismatch cannot replace valid symbols with NA.
target_anchored_reference_table <- read.csv(
  TARGET_ANCHORED_REFERENCE_FILE,
  stringsAsFactors = FALSE
)
if (!("Gene" %in% colnames(target_anchored_reference_table)) &&
    !("Raw_Gene" %in% colnames(target_anchored_reference_table))) {
  stop(
    "Target-anchored reference table lacks both Gene and Raw_Gene columns:\n  ",
    TARGET_ANCHORED_REFERENCE_FILE,
    call. = FALSE
  )
}
if (!("Raw_Gene" %in% colnames(target_anchored_reference_table))) {
  target_anchored_reference_table$Raw_Gene <-
    strip_M(target_anchored_reference_table$Gene)
}
target_anchored_reference_order <- unique(
  as.character(
    target_anchored_reference_table$Raw_Gene
  )
)
target_anchored_reference_order <-
  target_anchored_reference_order[
    !is.na(target_anchored_reference_order) &
      nzchar(target_anchored_reference_order)
  ]
canonical_target_signature <- add_M(
  target_anchored_reference_order
)
if (exists(
  "target_anchored_signature",
  envir = .GlobalEnv,
  inherits = FALSE
)) {
  loaded_signature <- get(
    "target_anchored_signature",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (!setequal(
    loaded_signature,
    canonical_target_signature
  )) {
    warning(
      "The target-anchored signature loaded by the common script did not ",
      "match the absolute discovery table. The canonical discovery-table ",
      "signature will be used for final-output plotting."
    )
  }
}
target_anchored_signature <- canonical_target_signature
message(">>> Project root:      ", PROJECT_ROOT)
message(">>> Results directory: ", RESULTS_DIR)
message(">>> Common script:     ", COMMON_SCRIPT_PATH)
message(">>> Input tables:      ", SENS_TAB_DIR)
message(">>> Final outputs:     ", SENS_FINAL_DIR)
message(
  ">>> Canonical target-anchored order: ",
  paste(
    target_anchored_reference_order,
    collapse = ", "
  )
)
has_patchwork <- requireNamespace(
  "patchwork",
  quietly = TRUE
)
has_ggrepel <- requireNamespace(
  "ggrepel",
  quietly = TRUE
)
has_patchwork_free <- has_patchwork &&
  "free" %in% getNamespaceExports("patchwork")
if (!requireNamespace(
  "scales",
  quietly = TRUE
)) {
  stop(
    "Package 'scales' is required for sensitivity figures.",
    call. = FALSE
  )
}
if (!has_patchwork) {
  stop(
    "Package 'patchwork' is required to assemble the main sensitivity figure.",
    call. = FALSE
  )
}
if (!has_ggrepel) {
  message(
    "Package 'ggrepel' is not installed; deterministic plain-text labels ",
    "will be used for SA6."
  )
}
if (!has_patchwork_free) {
  message(
    "patchwork::free() is unavailable; the script will use the compatible ",
    "fallback layout. Upgrade patchwork for the tightest panel-A alignment."
  )
}
# ------------------------------------------------------------------------------
# Paths and I/O helpers
# ------------------------------------------------------------------------------
sens_file <- function(filename) {
  file.path(
    SENS_TAB_DIR,
    filename
  )
}
private_file <- function(filename) {
  file.path(
    SENS_PRIVATE_DIR,
    filename
  )
}
final_file <- function(filename) {
  file.path(
    SENS_FINAL_DIR,
    filename
  )
}
read_sens <- function(filename, required = FALSE) {
  path <- sens_file(filename)
  if (!file.exists(path)) {
    if (required) {
      stop(
        "Required sensitivity output not found: ",
        path,
        call. = FALSE
      )
    }
    return(NULL)
  }
  tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      if (required) {
        stop(
          "Cannot read required sensitivity output: ",
          path,
          "\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
      warning(
        "Skipping unreadable sensitivity output: ",
        path
      )
      NULL
    }
  )
}
read_private <- function(filename, required = FALSE) {
  path <- private_file(filename)
  if (!file.exists(path)) {
    if (required) {
      stop(
        "Required participant-level sensitivity output not found: ",
        path,
        call. = FALSE
      )
    }
    return(NULL)
  }
  tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      if (required) {
        stop(
          "Cannot read required participant-level sensitivity output: ",
          path,
          "\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
      warning(
        "Skipping unreadable participant-level output: ",
        path
      )
      NULL
    }
  )
}
write_final <- function(x, filename) {
  out <- final_file(filename)
  write.csv(
    x,
    out,
    row.names = FALSE,
    na = ""
  )
  message(">>> Saved: ", out)
  invisible(out)
}
save_final_plot <- function(
    plot_obj,
    filename_base,
    width = 7,
    height = 5,
    dpi = 600
) {
  png_path <- final_file(
    paste0(filename_base, ".png")
  )
  pdf_path <- final_file(
    paste0(filename_base, ".pdf")
  )
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  message(">>> Saved: ", png_path)
  message(">>> Saved: ", pdf_path)
  invisible(
    c(
      png = png_path,
      pdf = pdf_path
    )
  )
}
has_cols <- function(df, cols) {
  !is.null(df) &&
    is.data.frame(df) &&
    nrow(df) > 0 &&
    all(cols %in% colnames(df))
}
safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}
safe_first <- function(df, col, default = NA) {
  if (is.null(df) ||
      !is.data.frame(df) ||
      nrow(df) == 0 ||
      !(col %in% colnames(df))) {
    return(default)
  }
  df[[col]][1]
}
safe_first_num <- function(
    df,
    col,
    default = NA_real_
) {
  suppressWarnings(
    as.numeric(
      safe_first(
        df,
        col,
        default = default
      )
    )
  )
}
safe_round <- function(x, digits = 3) {
  x <- suppressWarnings(
    as.numeric(x)
  )
  ifelse(
    is.finite(x),
    round(x, digits),
    NA_real_
  )
}
safe_sd_local2 <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x)
}
# ------------------------------------------------------------------------------
# Publication-oriented color palette
# ------------------------------------------------------------------------------
# The palette follows the earlier sensitivity figures: neutral greys for baseline
# or background, warm orange for target-anchored/primary results, and teal/blue
# for source/transfer diagnostics. It is intentionally restrained for BIB-style
# manuscript figures.
COL <- list(
  black = "#222222",
  dark_grey = "#555555",
  mid_grey = "#7A7A7A",
  light_grey = "#D9D9D9",
  pale_grey = "#F2F2F2",
  orange = "#D95F02",
  orange_light = "#F4A261",
  teal = "#1B9E9E",
  teal_light = "#80CDC1",
  blue = "#377EB8",
  blue_dark = "#2C5D85",
  purple = "#756BB1"
)
pal_metric <- c(
  "Target-anchored" = COL$orange,
  "Transfer-supported" = COL$teal,
  "Source-supported" = COL$blue_dark,
  "Core" = COL$purple
)
pal_dataset <- c(
  TARGET_NAME = COL$orange_light,
  SOURCE_NAME = COL$teal
)
names(pal_dataset) <- c(
  TARGET_NAME,
  SOURCE_NAME
)
class_palette <- c(
  "High-stability transfer core" = COL$orange,
  "Source-supported transfer-stabilized" = COL$teal,
  "Target-enhanced transfer-retained" = COL$blue,
  "Target-only source-discordant" = COL$purple,
  "Source-supported background" = COL$mid_grey,
  "Background" = COL$light_grey
)
class_shapes <- c(
  "High-stability transfer core" = 16,
  "Source-supported transfer-stabilized" = 17,
  "Target-enhanced transfer-retained" = 15,
  "Target-only source-discordant" = 3,
  "Source-supported background" = 16,
  "Background" = 16
)
# Intended full-width, two-column manuscript dimensions.
MAIN_FIG_WIDTH_IN <- 7.2
MAIN_FIG_HEIGHT_IN <- 6.6
MAIN_FIG_DPI <- 600
theme_vita <- function(base_size = 11) {
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        color = COL$black,
        size = base_size + 1,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        color = COL$dark_grey,
        size = base_size - 1,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 3)
      ),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(
        color = COL$black
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 4)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 4)
      ),
      axis.text = ggplot2::element_text(
        color = COL$black
      ),
      legend.title = ggplot2::element_text(
        color = COL$black
      ),
      legend.text = ggplot2::element_text(
        color = COL$black
      ),
      legend.margin = ggplot2::margin(
        0,
        0,
        0,
        0
      ),
      legend.box.spacing = grid::unit(
        2,
        "pt"
      ),
      strip.background = ggplot2::element_rect(
        fill = "white",
        color = COL$black
      ),
      strip.text = ggplot2::element_text(
        color = COL$black
      ),
      plot.margin = ggplot2::margin(
        5,
        5,
        5,
        5
      )
    )
}
# ------------------------------------------------------------------------------
# Read sensitivity outputs
# ------------------------------------------------------------------------------
dir.create(
  SENS_FINAL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)
sa1_sum <- read_sens(
  "SA1_NoiseKnockin_MainSummary.csv"
)
sa1_gene <- read_sens(
  "SA1_NoiseKnockin_GeneLevel.csv"
)
sa2_sum <- read_sens(
  "SA2_SourceDose_Summary.csv"
)
sa2_det <- read_sens(
  "SA2_SourceDose_Detailed.csv"
)
sa3_sum <- read_sens(
  "SA3_ThresholdGrid_Summary.csv"
)
sa3_gene <- read_sens(
  "SA3_ThresholdGrid_GeneLevel.csv"
)
sa4_sum <- read_sens(
  "SA4_FixedSignature_Prediction_Summary.csv"
)
sa4_all <- read_sens(
  "SA4_FixedSignature_RepeatedCV_All.csv"
)
sa4_contrast <- read_sens(
  "SA4_FixedSignature_Pairwise_Contrasts.csv"
)
sa5_sum <- read_sens(
  "SA5_CoxSnell_Residual_Summary.csv"
)
sa5_data <- read_sens(
  "SA5_CoxSnell_Residual_Data.csv"
)
sa6_weights <- read_sens(
  "SA6_SourceWeights_FullData.csv"
)
sa6_pca_raw <- read_private(
  "SA6_PCA_SourceTarget_RawPooled_Data.csv"
)
if (is.null(sa6_pca_raw)) {
  sa6_pca_raw <- read_private(
    "SA6_PCA_SourceTarget_Data.csv"
  )
}
sa6_pca_within <- read_private(
  "SA6_PCA_SourceTarget_WithinCohortZ_Data.csv"
)
sa6_stab <- read_sens(
  "SA6_SourcePIP_TargetStability_Data.csv"
)
sa6_corr <- read_sens(
  "SA6_SourcePIP_TargetStability_Correlation.csv"
)
sa6_quad <- read_sens(
  "SA6_SourceTarget_Quadrant_Summary.csv"
)
sa7_sum <- read_sens(
  "SA7_HyperparameterGrid_Summary.csv"
)
sa7_gene <- read_sens(
  "SA7_HyperparameterGrid_GeneLevel.csv"
)
sa8_rule <- read_sens(
  "SA8_TargetAnchoring_RuleSummary.csv"
)
sa8_gene <- read_sens(
  "SA8_TargetAnchoring_GeneMembership.csv"
)
sa8_cv <- read_sens(
  "SA8_TargetAnchoring_CIndexComparison.csv"
)
sa8_eff <- read_sens(
  "SA8_TargetAnchoring_EffectComparison.csv"
)
# Main-text panels must be complete. Supplementary inputs remain optional so
# that a missing diagnostic does not invalidate already completed main analyses.
main_inputs <- list(
  SA1_main_summary = sa1_sum,
  SA1_gene_level = sa1_gene,
  SA2_source_dose_summary = sa2_sum,
  SA3_threshold_grid_summary = sa3_sum,
  SA4_fixed_signature_summary = sa4_sum
)
missing_main_inputs <- names(main_inputs)[vapply(
  main_inputs,
  function(x) {
    is.null(x) ||
      !is.data.frame(x) ||
      nrow(x) == 0
  },
  logical(1)
)]
if (length(missing_main_inputs) > 0) {
  stop(
    "Cannot generate the complete main sensitivity figure. Missing inputs: ",
    paste(
      missing_main_inputs,
      collapse = ", "
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------------------------
# Final table outputs
# ------------------------------------------------------------------------------
if (has_cols(
  sa1_gene,
  c(
    "In_TargetAnchored_Reference",
    "In_TransferSupported_Reference",
    "In_SourceSupported_Reference",
    "In_Core_Reference",
    "In_TargetOnly_Reference",
    "In_TargetDiscordant_Reference",
    "Is_Fake",
    "Selected_TargetAnchored",
    "VITA_Stability"
  )
)) {
  sa1_table <- sa1_gene %>%
    dplyr::filter(
      In_TargetAnchored_Reference |
        In_TransferSupported_Reference |
        In_SourceSupported_Reference |
        In_Core_Reference |
        In_TargetOnly_Reference |
        In_TargetDiscordant_Reference |
        Is_Fake |
        Selected_TargetAnchored
    ) %>%
    dplyr::arrange(
      dplyr::desc(Selected_TargetAnchored),
      dplyr::desc(Is_Fake),
      dplyr::desc(VITA_Stability)
    )
  write_final(
    sa1_table,
    "Table_SA1_NoiseKnockin_Final.csv"
  )
} else if (!is.null(sa1_gene)) {
  write_final(
    sa1_gene,
    "Table_SA1_NoiseKnockin_Final.csv"
  )
}
if (!is.null(sa2_sum)) {
  write_final(
    sa2_sum,
    "Table_SA2_SourceDose_Final.csv"
  )
}
if (!is.null(sa3_sum)) {
  write_final(
    sa3_sum,
    "Table_SA3_ThresholdGrid_Final.csv"
  )
}
if (!is.null(sa3_gene)) {
  write_final(
    sa3_gene,
    "Table_SA3_ThresholdGrid_GeneLevel_Final.csv"
  )
}
if (!is.null(sa4_sum)) {
  write_final(
    sa4_sum,
    "Table_SA4_FixedSignature_Prediction_Final.csv"
  )
}
if (!is.null(sa4_contrast)) {
  write_final(
    sa4_contrast,
    "Table_SA4_FixedSignature_Contrasts_Final.csv"
  )
}
if (!is.null(sa5_sum)) {
  write_final(
    sa5_sum,
    "Table_SA5_CoxSnell_Final.csv"
  )
}
if (!is.null(sa6_weights)) {
  write_final(
    sa6_weights,
    "Table_SA6_SourceWeights_Final.csv"
  )
}
# Participant-level PCA coordinate tables are used to draw SA6 figures but are
# deliberately not copied into final_outputs.
if (!is.null(sa6_stab)) {
  write_final(
    sa6_stab,
    "Table_SA6_SourcePIP_TargetStability_Final.csv"
  )
}
if (!is.null(sa6_corr)) {
  write_final(
    sa6_corr,
    "Table_SA6_SourcePIP_TargetStability_Correlation_Final.csv"
  )
}
if (!is.null(sa6_quad)) {
  write_final(
    sa6_quad,
    "Table_SA6_SourceTarget_Quadrant_Final.csv"
  )
}
if (!is.null(sa7_sum)) {
  write_final(
    sa7_sum,
    "Table_SA7_HyperparameterGrid_Final.csv"
  )
}
if (!is.null(sa7_gene)) {
  write_final(
    sa7_gene,
    "Table_SA7_HyperparameterGrid_GeneLevel_Final.csv"
  )
}
if (!is.null(sa8_rule)) {
  write_final(
    sa8_rule,
    "Table_SA8_TargetAnchoringRuleAudit_Final.csv"
  )
}
if (!is.null(sa8_gene)) {
  write_final(
    sa8_gene,
    "Table_SA8_TargetAnchoringGeneMembership_Final.csv"
  )
}
if (!is.null(sa8_cv)) {
  write_final(
    sa8_cv,
    "Table_SA8_TargetAnchoringCIndexComparison_Final.csv"
  )
}
if (!is.null(sa8_eff)) {
  write_final(
    sa8_eff,
    "Table_SA8_TargetAnchoringEffectComparison_Final.csv"
  )
}
# ------------------------------------------------------------------------------
# Unified sensitivity summary table
# ------------------------------------------------------------------------------
summary_rows <- list()
if (!is.null(sa1_sum) &&
    nrow(sa1_sum) > 0) {
  summary_rows[["SA1"]] <- data.frame(
    Analysis = "SA1 Noise knock-in",
    Main_Result = paste0(
      "Synthetic noise selected=",
      safe_first(
        sa1_sum,
        "Noise_Selected"
      ),
      "; noise selection rate=",
      safe_round(
        safe_first_num(
          sa1_sum,
          "Noise_Selection_Rate"
        ),
        3
      ),
      "; target-anchored retention=",
      safe_round(
        safe_first_num(
          sa1_sum,
          "TargetAnchored_Recall"
        ),
        3
      ),
      "; transfer-supported retention=",
      safe_round(
        safe_first_num(
          sa1_sum,
          "TransferSupported_Recall"
        ),
        3
      ),
      "; target-discordant retention=",
      safe_round(
        safe_first_num(
          sa1_sum,
          "TargetDiscordant_Retention"
        ),
        3
      )
    ),
    Interpretation = paste(
      "Noise knock-in evaluates whether the target-anchored",
      "reported set and its subclasses remain stable after",
      "adding synthetic null features."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa2_sum,
  c(
    "Dose",
    "Avg_C_Index",
    "Avg_TargetAnchored_Recall",
    "Avg_TransferSupported_Recall"
  )
)) {
  summary_rows[["SA2"]] <- data.frame(
    Analysis = "SA2 Source dose-response",
    Main_Result = paste(
      paste0(
        sa2_sum$Dose,
        ": C=",
        safe_round(
          sa2_sum$Avg_C_Index,
          3
        ),
        ", target-anchored recall=",
        safe_round(
          sa2_sum$Avg_TargetAnchored_Recall,
          3
        ),
        ", transfer-supported recall=",
        safe_round(
          sa2_sum$Avg_TransferSupported_Recall,
          3
        ),
        if ("Avg_TargetDiscordant_Retention" %in%
            colnames(sa2_sum)) {
          paste0(
            ", target-discordant retention=",
            safe_round(
              sa2_sum$Avg_TargetDiscordant_Retention,
              3
            )
          )
        } else {
          ""
        }
      ),
      collapse = "; "
    ),
    Interpretation = paste(
      "Source dose-response uses repeated downsampling at",
      "30% and 60%, with 100% shown as the full-source reference."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa3_sum,
  c(
    "TARGET_BFDR",
    "STABILITY_CUTOFF",
    "N_TargetAnchored_Selected",
    "TargetAnchored_Recall",
    "TransferSupported_Recall",
    "TargetDiscordant_Retention"
  )
)) {
  main_cell <- sa3_sum %>%
    dplyr::filter(
      abs(
        safe_num(TARGET_BFDR) -
          get(
            "TARGET_BFDR",
            envir = .GlobalEnv
          )
      ) < 1e-8,
      abs(
        safe_num(STABILITY_CUTOFF) -
          get(
            "EXTENDED_STABILITY_CUTOFF",
            envir = .GlobalEnv
          )
      ) < 1e-8
    )
  if (nrow(main_cell) > 0) {
    main_result <- paste0(
      "Main grid cell selected ",
      main_cell$N_TargetAnchored_Selected[1],
      " target-anchored genes; target-anchored recall=",
      safe_round(
        main_cell$TargetAnchored_Recall[1],
        3
      ),
      "; transfer-supported recall=",
      safe_round(
        main_cell$TransferSupported_Recall[1],
        3
      ),
      "; target-discordant retention=",
      safe_round(
        main_cell$TargetDiscordant_Retention[1],
        3
      )
    )
  } else {
    main_result <- paste(
      "Threshold grid completed; the prespecified",
      "main grid cell was not found."
    )
  }
  summary_rows[["SA3"]] <- data.frame(
    Analysis = "SA3 BFDR/stability threshold grid",
    Main_Result = main_result,
    Interpretation = paste(
      "The grid characterizes whether the target-anchored",
      "reported set depends strongly on BFDR and",
      "bootstrap-stability thresholds."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa4_sum,
  c(
    "Model",
    "RepeatedCV_Mean_C_Index",
    "Mean_Delta_vs_Clinical",
    "Positive_Delta_Rate"
  )
)) {
  rows_sig <- sa4_sum %>%
    dplyr::filter(
      Model != "Clinical-Only"
    )
  summary_rows[["SA4"]] <- data.frame(
    Analysis = "SA4 Fixed-signature prediction",
    Main_Result = paste(
      paste0(
        rows_sig$Model,
        ": C=",
        safe_round(
          rows_sig$RepeatedCV_Mean_C_Index,
          3
        ),
        ", Delta C=",
        safe_round(
          rows_sig$Mean_Delta_vs_Clinical,
          3
        ),
        ", positive rate=",
        safe_round(
          rows_sig$Positive_Delta_Rate,
          3
        )
      ),
      collapse = "; "
    ),
    Interpretation = paste(
      "Fixed-signature internal repeated CV compares",
      "target-anchored, transfer-supported, target-only,",
      "and source-discordant gene sets."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa5_sum,
  c(
    "Signature_Set",
    "RMSE_from_UnitSlope",
    "MaxAbs_Deviation"
  )
)) {
  summary_rows[["SA5"]] <- data.frame(
    Analysis = "SA5 Cox-Snell AFT diagnostic",
    Main_Result = paste(
      paste0(
        sa5_sum$Signature_Set,
        ": RMSE=",
        safe_round(
          sa5_sum$RMSE_from_UnitSlope,
          3
        ),
        ", max deviation=",
        safe_round(
          sa5_sum$MaxAbs_Deviation,
          3
        )
      ),
      collapse = "; "
    ),
    Interpretation = paste(
      "Residual behavior is used as a log-normal AFT",
      "assumption diagnostic for fixed reported gene sets."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa6_corr,
  c(
    "Metric",
    "Estimate"
  )
)) {
  summary_rows[["SA6"]] <- data.frame(
    Analysis = "SA6 Source-target transferability",
    Main_Result = paste(
      paste0(
        sa6_corr$Metric,
        "=",
        safe_round(
          sa6_corr$Estimate,
          3
        )
      ),
      collapse = "; "
    ),
    Interpretation = paste(
      "Transferability diagnostics evaluate raw/standardized",
      "source-target feature distributions and source PIP",
      "versus target-domain stability."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa7_sum,
  c(
    "TargetAnchored_Recall",
    "TargetDiscordant_Retention"
  )
)) {
  successful <- if ("Fit_Success" %in%
                    colnames(sa7_sum)) {
    sum(
      sa7_sum$Fit_Success %in% TRUE,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }
  summary_rows[["SA7"]] <- data.frame(
    Analysis = "SA7 Hyperparameter grid",
    Main_Result = paste0(
      "Completed ",
      nrow(sa7_sum),
      " reduced-grid parameter cells; successful cells=",
      successful,
      "; median target-anchored recall=",
      safe_round(
        stats::median(
          safe_num(
            sa7_sum$TargetAnchored_Recall
          ),
          na.rm = TRUE
        ),
        3
      ),
      "; median target-discordant retention=",
      safe_round(
        stats::median(
          safe_num(
            sa7_sum$TargetDiscordant_Retention
          ),
          na.rm = TRUE
        ),
        3
      )
    ),
    Interpretation = paste(
      "The reduced grid evaluates whether target-anchored",
      "and transfer-supported discoveries depend strongly",
      "on prior sparsity, transfer strength, tau, or",
      "beta-strength transfer."
    ),
    stringsAsFactors = FALSE
  )
}
if (has_cols(
  sa8_rule,
  c(
    "Rule_Set",
    "N_Genes"
  )
)) {
  summary_rows[["SA8"]] <- data.frame(
    Analysis = "SA8 Target-anchoring rule audit",
    Main_Result = paste(
      paste0(
        sa8_rule$Rule_Set,
        ": n=",
        sa8_rule$N_Genes
      ),
      collapse = "; "
    ),
    Interpretation = paste(
      "This audit separates transfer-supported genes from",
      "target-only source-discordant genes, preventing",
      "source-domain evidence from acting as a veto."
    ),
    stringsAsFactors = FALSE
  )
}
summary_table <- if (length(summary_rows) > 0) {
  dplyr::bind_rows(summary_rows)
} else {
  data.frame(
    Analysis = character(0),
    Main_Result = character(0),
    Interpretation = character(0),
    stringsAsFactors = FALSE
  )
}
write_final(
  summary_table,
  "Table_Glioma_Sensitivity_Overall_Summary.csv"
)
print(summary_table)
# ------------------------------------------------------------------------------
# Helper functions for plotting
# ------------------------------------------------------------------------------
model_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "VITA_HDM_TargetAnchored_Reported_Set" =
      "Target-anchored",
    "TargetAnchored_Reported_Set" =
      "Target-anchored",
    "VITA_HDM_TransferSupported_Set" =
      "Transfer-supported",
    "TransferSupported_Set" =
      "Transfer-supported",
    "VITA_HDM_SourceSupported_Subset" =
      "Source-supported subset",
    "VITA_HDM_Core" =
      "High-stability core",
    "TargetOnly_Stable" =
      "Target-only stable",
    "TargetOnly_SourceDiscordant" =
      "Target-only source-discordant",
    "Clinical-Only" =
      "Clinical-only",
    .default = as.character(x)
  )
}
# ------------------------------------------------------------------------------
# Panel A: noise knock-in robustness table
# ------------------------------------------------------------------------------
make_panel_A_noise <- function() {
  required_cols <- c(
    "Gene",
    "Raw_Gene",
    "In_TargetAnchored_Reference",
    "Selected_TargetAnchored",
    "VITA_Stability",
    "TargetOnly_Stability"
  )
  if (!has_cols(
    sa1_gene,
    required_cols
  )) {
    return(NULL)
  }
  baseline_df <- read_sens(
    "Baseline_TargetAnchored_Classification.csv",
    required = TRUE
  )
  baseline_required <- c(
    "Gene",
    "VITA_Stability",
    "TargetOnly_Stability",
    "Source_Supported",
    "TargetOnly_Stable",
    "TargetOnly_SourceDiscordant",
    "Target_Enhanced",
    "Category"
  )
  if (!has_cols(
    baseline_df,
    baseline_required
  )) {
    stop(
      "Baseline_TargetAnchored_Classification.csv lacks columns required ",
      "for Figure 4A.",
      call. = FALSE
    )
  }
  normalize_gene_key <- function(x) {
    x <- as.character(x)
    x <- ifelse(
      is.na(x),
      "",
      x
    )
    add_M(
      strip_M(x)
    )
  }
  baseline_join <- baseline_df %>%
    dplyr::transmute(
      Gene_Key = normalize_gene_key(Gene),
      Baseline_VITA_Stability =
        safe_num(VITA_Stability),
      Baseline_TargetOnly_Stability =
        safe_num(TargetOnly_Stability),
      Baseline_Source_Supported =
        Source_Supported %in% TRUE,
      Baseline_TargetOnly_Stable =
        TargetOnly_Stable %in% TRUE,
      Baseline_TargetOnly_SourceDiscordant =
        TargetOnly_SourceDiscordant %in% TRUE,
      Baseline_Target_Enhanced =
        Target_Enhanced %in% TRUE,
      Baseline_Category =
        as.character(Category)
    ) %>%
    dplyr::distinct(
      Gene_Key,
      .keep_all = TRUE
    )
  tab <- sa1_gene %>%
    dplyr::filter(
      In_TargetAnchored_Reference %in% TRUE
    ) %>%
    dplyr::mutate(
      Gene_Key = normalize_gene_key(Gene),
      Raw_Gene = dplyr::if_else(
        is.na(Raw_Gene) |
          !nzchar(as.character(Raw_Gene)),
        strip_M(Gene_Key),
        as.character(Raw_Gene)
      )
    ) %>%
    dplyr::left_join(
      baseline_join,
      by = "Gene_Key"
    )
  if (nrow(tab) == 0) {
    return(NULL)
  }
  missing_baseline <- tab$Raw_Gene[
    is.na(tab$Baseline_VITA_Stability) &
      is.na(tab$Baseline_TargetOnly_Stability)
  ]
  if (length(missing_baseline) > 0) {
    stop(
      "Figure 4A baseline matching failed for: ",
      paste(
        unique(missing_baseline),
        collapse = ", "
      ),
      "\nCheck that SA1 and baseline tables came from the same results directory.",
      call. = FALSE
    )
  }
  # Use target-only stability for target-enhanced/source-discordant genes and
  # VITA-HDM stability for source-supported transfer genes.
  tab <- tab %>%
    dplyr::mutate(
      Use_TargetOnly_Score =
        Baseline_TargetOnly_SourceDiscordant |
        Baseline_Target_Enhanced,
      Baseline_Score = dplyr::if_else(
        Use_TargetOnly_Score,
        Baseline_TargetOnly_Stability,
        Baseline_VITA_Stability
      ),
      Perturbed_Score = dplyr::if_else(
        Use_TargetOnly_Score,
        safe_num(TargetOnly_Stability),
        safe_num(VITA_Stability)
      ),
      Status = dplyr::if_else(
        Selected_TargetAnchored %in% TRUE,
        "Retained",
        "Not retained"
      ),
      Baseline = ifelse(
        is.finite(Baseline_Score),
        sprintf(
          "%.2f",
          Baseline_Score
        ),
        "--"
      ),
      Perturbed = ifelse(
        is.finite(Perturbed_Score),
        sprintf(
          "%.2f",
          Perturbed_Score
        ),
        "--"
      )
    )
  # Keep gene symbols as character values; use a separate numeric display order.
  display_order <- unique(c(
    intersect(
      target_anchored_reference_order,
      tab$Raw_Gene
    ),
    setdiff(
      tab$Raw_Gene,
      target_anchored_reference_order
    )
  ))
  tab <- tab %>%
    dplyr::mutate(
      Display_Order = match(
        Raw_Gene,
        display_order
      )
    ) %>%
    dplyr::arrange(
      Display_Order
    )
  if (anyNA(tab$Raw_Gene) ||
      any(!nzchar(tab$Raw_Gene))) {
    stop(
      "Figure 4A contains an empty or NA gene symbol after ordering.",
      call. = FALSE
    )
  }
  if (anyDuplicated(tab$Raw_Gene)) {
    stop(
      "Figure 4A contains duplicated reference gene symbols: ",
      paste(
        unique(
          tab$Raw_Gene[
            duplicated(tab$Raw_Gene)
          ]
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  table_out <- tab %>%
    dplyr::transmute(
      Gene = Raw_Gene,
      Baseline = Baseline,
      Perturbed = Perturbed,
      Status = Status,
      Baseline_Category = Baseline_Category
    )
  write_final(
    table_out,
    "Table_Figure4A_NoiseKnockin_MainText.csv"
  )
  n <- nrow(tab)
  # Keep the table left-aligned with the panel title while reducing its visual
  # footprint. The extra space on the right and below prevents Panel A from
  # dominating the upper row of the composite figure.
  table_col_centers <- c(
    0.50,
    1.40,
    2.30,
    3.20
  )
  table_col_half_width <- 0.45
  table_row_step <- 0.84
  table_row_half_height <- 0.40
  table_header_y <- (n + 1) * table_row_step
  tab$Row <- rev(
    seq_len(n)
  ) * table_row_step
  header <- data.frame(
    Col = c(
      "Gene",
      "Baseline",
      "Perturbed",
      "Status"
    ),
    X = table_col_centers,
    stringsAsFactors = FALSE
  )
  row_rects <- expand.grid(
    Row_Index = seq_len(n),
    X = header$X
  )
  row_rects$Row <-
    tab$Row[row_rects$Row_Index]
  row_rects$Fill <- ifelse(
    row_rects$Row_Index %% 2 == 0,
    COL$pale_grey,
    "white"
  )
  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = row_rects,
      ggplot2::aes(
        xmin = X - table_col_half_width,
        xmax = X + table_col_half_width,
        ymin = Row - table_row_half_height,
        ymax = Row + table_row_half_height,
        fill = Fill
      ),
      color = "white",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_rect(
      data = header,
      ggplot2::aes(
        xmin = X - table_col_half_width,
        xmax = X + table_col_half_width,
        ymin = table_header_y -
          table_row_half_height,
        ymax = table_header_y +
          table_row_half_height
      ),
      fill = COL$dark_grey,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      data = header,
      ggplot2::aes(
        x = X,
        y = table_header_y,
        label = Col
      ),
      color = "white",
      fontface = "bold",
      size = 3.10
    ) +
    ggplot2::geom_text(
      data = tab,
      ggplot2::aes(
        y = Row,
        label = Raw_Gene
      ),
      x = table_col_centers[1],
      fontface = "italic",
      size = 3.00
    ) +
    ggplot2::geom_text(
      data = tab,
      ggplot2::aes(
        y = Row,
        label = Baseline
      ),
      x = table_col_centers[2],
      size = 3.00
    ) +
    ggplot2::geom_text(
      data = tab,
      ggplot2::aes(
        y = Row,
        label = Perturbed
      ),
      x = table_col_centers[3],
      size = 3.00
    ) +
    ggplot2::geom_text(
      data = tab,
      ggplot2::aes(
        y = Row,
        label = Status,
        color = Status
      ),
      x = table_col_centers[4],
      size = 3.00,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Retained" = COL$black,
        "Not retained" = COL$orange
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0.05,
        4.15
      ),
      ylim = c(
        -0.45,
        table_header_y +
          table_row_half_height
      ),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "A. Noise knock-in robustness",
      subtitle = paste0(
        "0/",
        safe_first(
          sa1_sum,
          "Injected_Noise_Features",
          "50"
        ),
        " noise genes selected; ",
        round(
          100 *
            safe_first_num(
              sa1_sum,
              "TargetAnchored_Recall"
            )
        ),
        "% target-anchored retention"
      )
    ) +
    ggplot2::theme_void(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        color = COL$black,
        size = 11,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 2)
      ),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(
        color = COL$dark_grey,
        size = 9,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 3)
      ),
      plot.margin = ggplot2::margin(
        3,
        3,
        3,
        3
      )
    )
}
# ------------------------------------------------------------------------------
# Panel B: source dose-response recovery
# ------------------------------------------------------------------------------
make_panel_B_source_recovery <- function() {
  required_cols <- c(
    "Dose",
    "Avg_TargetAnchored_Recall",
    "Avg_TransferSupported_Recall",
    "Avg_SourceSupported_Recall"
  )
  if (!has_cols(
    sa2_sum,
    required_cols
  )) {
    return(NULL)
  }
  df <- sa2_sum
  df$Dose_Num <- safe_num(
    gsub(
      "%",
      "",
      df$Dose
    )
  )
  if (!("Analysis_Type" %in% colnames(df))) {
    df$Analysis_Type <- ifelse(
      df$Dose_Num == 100,
      "Full-source reference",
      "Reduced-source repeated subsampling"
    )
  }
  df <- df %>%
    dplyr::filter(
      is.finite(Dose_Num)
    ) %>%
    dplyr::arrange(
      Dose_Num
    )
  recovery_map <- c(
    "Source-supported" =
      "Avg_SourceSupported_Recall",
    "Transfer-supported" =
      "Avg_TransferSupported_Recall",
    "Target-anchored" =
      "Avg_TargetAnchored_Recall"
  )
  recovery_map <- recovery_map[
    recovery_map %in% colnames(df)
  ]
  if (length(recovery_map) == 0) {
    return(NULL)
  }
  bar_df <- dplyr::bind_rows(
    lapply(
      names(recovery_map),
      function(metric) {
        col <- recovery_map[[metric]]
        sd_col <- paste0(
          "SD_",
          sub(
            "^Avg_",
            "",
            col
          )
        )
        data.frame(
          Dose_Num = df$Dose_Num,
          Dose_Label = paste0(
            df$Dose_Num,
            "%"
          ),
          Metric = metric,
          Rate = safe_num(
            df[[col]]
          ),
          SD = if (sd_col %in%
                   colnames(df)) {
            safe_num(
              df[[sd_col]]
            )
          } else {
            NA_real_
          },
          Analysis_Type = df$Analysis_Type,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  dose_levels <- paste0(
    sort(
      unique(
        bar_df$Dose_Num
      )
    ),
    "%"
  )
  bar_df <- bar_df %>%
    dplyr::mutate(
      Dose_Label = factor(
        Dose_Label,
        levels = dose_levels
      ),
      Dose_Index = as.numeric(
        Dose_Label
      ),
      Metric = factor(
        Metric,
        levels = c(
          "Source-supported",
          "Transfer-supported",
          "Target-anchored"
        )
      ),
      Is_Reference =
        grepl(
          "Full-source",
          Analysis_Type
        ) |
        Dose_Num == 100,
      ymin = pmax(
        0,
        Rate -
          ifelse(
            is.finite(SD),
            SD,
            0
          )
      ),
      ymax = pmin(
        1,
        Rate +
          ifelse(
            is.finite(SD),
            SD,
            0
          )
      )
    )
  ggplot2::ggplot(
    bar_df,
    ggplot2::aes(
      x = Dose_Index,
      y = Rate,
      fill = Metric
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(
        width = 0.72
      ),
      width = 0.62,
      alpha = 0.88,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_errorbar(
      data = bar_df %>%
        dplyr::filter(
          !Is_Reference,
          is.finite(SD)
        ),
      ggplot2::aes(
        ymin = ymin,
        ymax = ymax,
        color = Metric
      ),
      position = ggplot2::position_dodge(
        width = 0.72
      ),
      width = 0.18,
      linewidth = 0.55,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = pal_metric,
      name = NULL
    ) +
    ggplot2::scale_color_manual(
      values = pal_metric,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(
        dose_levels
      ),
      labels = dose_levels,
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        1.04
      ),
      labels = scales::label_percent(
        accuracy = 1
      ),
      breaks = seq(
        0,
        1,
        0.25
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0.45,
        length(dose_levels) + 0.45
      ),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    ) +
    ggplot2::labs(
      title = "B. Source dose-response recovery",
      subtitle = paste(
        "30/60% subsampling;",
        "100% full-source reference"
      ),
      x = "Source data proportion",
      y = "Gene-set recovery rate"
    ) +
    theme_vita(
      base_size = 10
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.box = "horizontal",
      legend.box.spacing = grid::unit(
        1,
        "pt"
      ),
      legend.margin = ggplot2::margin(
        t = 1
      ),
      legend.key.width = grid::unit(
        0.65,
        "lines"
      ),
      legend.key.height = grid::unit(
        0.65,
        "lines"
      ),
      legend.spacing.x = grid::unit(
        0.12,
        "lines"
      ),
      legend.spacing.y = grid::unit(
        0.05,
        "lines"
      ),
      legend.text = ggplot2::element_text(
        size = 7.8
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 3)
      ),
      axis.text.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 2)
      ),
      plot.margin = ggplot2::margin(
        3,
        5,
        3,
        3
      )
    )
}
# ------------------------------------------------------------------------------
# Panel D: repeated-CV gain of fixed signatures
# ------------------------------------------------------------------------------
make_panel_D_prediction_gain <- function() {
  if (!has_cols(
    sa4_sum,
    c(
      "Model",
      "Mean_Delta_vs_Clinical"
    )
  )) {
    return(NULL)
  }
  df <- sa4_sum %>%
    dplyr::filter(
      Model != "Clinical-Only"
    ) %>%
    dplyr::mutate(
      Model_Label_Character =
        model_label(Model),
      Delta =
        safe_num(Mean_Delta_vs_Clinical)
    ) %>%
    dplyr::filter(
      is.finite(Delta)
    )
  if (nrow(df) == 0) {
    return(NULL)
  }
  preferred <- c(
    "High-stability core",
    "Source-supported subset",
    "Target-only source-discordant",
    "Target-anchored",
    "Transfer-supported",
    "Target-only stable"
  )
  if (anyNA(df$Model_Label_Character) ||
      any(!nzchar(
        df$Model_Label_Character
      ))) {
    stop(
      "Panel D contains an empty or NA publication label.",
      call. = FALSE
    )
  }
  if (anyDuplicated(
    df$Model_Label_Character
  )) {
    stop(
      "Panel D contains duplicated publication labels: ",
      paste(
        unique(
          df$Model_Label_Character[
            duplicated(
              df$Model_Label_Character
            )
          ]
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  observed_labels <- unique(
    df$Model_Label_Character
  )
  unexpected_labels <- setdiff(
    observed_labels,
    preferred
  )
  if (length(unexpected_labels) > 0) {
    warning(
      "Panel D encountered unlisted model labels; they will be appended: ",
      paste(
        unexpected_labels,
        collapse = ", "
      )
    )
  }
  plot_levels <- c(
    intersect(
      preferred,
      observed_labels
    ),
    unexpected_labels
  )
  df <- df %>%
    dplyr::mutate(
      Model_Label = factor(
        Model_Label_Character,
        levels = plot_levels
      )
    )
  if (anyNA(df$Model_Label)) {
    stop(
      "Panel D factor conversion produced an NA label.",
      call. = FALSE
    )
  }
  # The top three bars are highlighted. All six labels use the same font size
  # and weight; only their color and placement differ.
  highlighted_labels <- c(
    "Target-anchored",
    "Transfer-supported",
    "Target-only stable"
  )
  positive_max <- max(
    c(
      0,
      df$Delta
    ),
    na.rm = TRUE
  )
  
  label_offset <- max(
    0.0006,
    0.015 * positive_max
  )
  bar_label_size <- 2.80
  bar_label_face <- "plain"
  df <- df %>%
    dplyr::mutate(
      Label_Placement = dplyr::if_else(
        Model_Label_Character %in%
          highlighted_labels,
        "inside",
        "outside"
      ),
      Label_Position = dplyr::if_else(
        Label_Placement == "inside",
        label_offset,
        Delta + label_offset
      )
    )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Model_Label,
      y = Delta
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = COL$light_grey,
      linewidth = 0.45
    ) +
    ggplot2::geom_col(
      width = 0.66,
      fill = COL$mid_grey
    ) +
    ggplot2::geom_col(
      data = df %>%
        dplyr::filter(
          Model_Label_Character %in%
            highlighted_labels
        ),
      width = 0.66,
      ggplot2::aes(
        fill = Model_Label_Character
      ),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = df %>%
        dplyr::filter(
          Label_Placement == "inside"
        ),
      ggplot2::aes(
        y = Label_Position,
        label = Model_Label_Character
      ),
      hjust = 0,
      color = "white",
      fontface = bar_label_face,
      size = bar_label_size,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = df %>%
        dplyr::filter(
          Label_Placement == "outside"
        ),
      ggplot2::aes(
        y = Label_Position,
        label = Model_Label_Character
      ),
      hjust = 0,
      color = COL$black,
      fontface = bar_label_face,
      size = bar_label_size,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Target-anchored" = COL$orange,
        "Transfer-supported" = COL$teal,
        "Target-only stable" = COL$blue
      )
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(
        add = 0.35
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.01,
          0.12
        )
      )
    ) +
    ggplot2::coord_flip(
      clip = "off"
    ) +
    ggplot2::labs(
      title = "D. Repeated-CV gain of fixed signatures",
      x = NULL,
      y = "Mean \u0394 C-index vs clinical-only"
    ) +
    theme_vita(
      base_size = 10
    ) +
    ggplot2::theme(
      axis.text.y =
        ggplot2::element_blank(),
      axis.ticks.y =
        ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        3,
        5,
        3,
        3
      )
    )
}
# ------------------------------------------------------------------------------
# Supplementary source-support versus target-stability plot
# ------------------------------------------------------------------------------
make_panel_D_source_pip <- function(
    include_legend = FALSE
) {
  if (!has_cols(
    sa6_stab,
    c(
      "TargetOnly_Stability",
      "Source_PIP_Weighted",
      "TargetAnchored_Reported",
      "Raw_Gene"
    )
  )) {
    return(NULL)
  }
  class_label_short <- c(
    "High-stability transfer core" =
      "Core",
    "Source-supported transfer-stabilized" =
      "Source-supported",
    "Target-enhanced transfer-retained" =
      "Target-enhanced",
    "Target-only source-discordant" =
      "Source-discordant",
    "Source-supported background" =
      "Source-supported background",
    "Background" =
      "Background"
  )
  class_palette_short <- c(
    "Core" = COL$orange,
    "Source-supported" = COL$teal,
    "Target-enhanced" = COL$blue,
    "Source-discordant" = COL$purple,
    "Source-supported background" =
      COL$mid_grey,
    "Background" = COL$light_grey
  )
  class_shapes_short <- c(
    "Core" = 16,
    "Source-supported" = 17,
    "Target-enhanced" = 15,
    "Source-discordant" = 3,
    "Source-supported background" = 16,
    "Background" = 16
  )
  df <- sa6_stab %>%
    dplyr::mutate(
      TargetOnly_Stability =
        safe_num(TargetOnly_Stability),
      Source_PIP_Weighted =
        safe_num(Source_PIP_Weighted),
      Plot_Class = ifelse(
        "Plot_Class" %in% colnames(.),
        as.character(Plot_Class),
        "Background"
      ),
      Plot_Class_Short = dplyr::recode(
        Plot_Class,
        !!!class_label_short,
        .default = Plot_Class
      )
    ) %>%
    dplyr::filter(
      is.finite(TargetOnly_Stability),
      is.finite(Source_PIP_Weighted)
    )
  if (nrow(df) == 0) {
    return(NULL)
  }
  bg <- df %>%
    dplyr::filter(
      !TargetAnchored_Reported
    )
  lab <- df %>%
    dplyr::filter(
      TargetAnchored_Reported
    )
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = SOURCE_SUPPORT_CUTOFF,
      linetype = "dashed",
      linewidth = 0.7,
      color = COL$black
    ) +
    ggplot2::geom_vline(
      xintercept =
        TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF,
      linetype = "dashed",
      linewidth = 0.7,
      color = COL$black
    ) +
    ggplot2::geom_point(
      data = bg,
      ggplot2::aes(
        x = TargetOnly_Stability,
        y = Source_PIP_Weighted
      ),
      color = COL$mid_grey,
      alpha = 0.18,
      size = 1.1
    ) +
    ggplot2::geom_point(
      data = lab,
      ggplot2::aes(
        x = TargetOnly_Stability,
        y = Source_PIP_Weighted,
        color = Plot_Class_Short,
        shape = Plot_Class_Short
      ),
      size = 2.7,
      alpha = 0.98
    ) +
    ggplot2::scale_color_manual(
      values = class_palette_short,
      na.value = COL$black,
      breaks = c(
        "Core",
        "Source-supported",
        "Target-enhanced",
        "Source-discordant"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = class_shapes_short,
      na.value = 16,
      breaks = c(
        "Core",
        "Source-supported",
        "Target-enhanced",
        "Source-discordant"
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 3
        )
      ),
      shape = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 3
        )
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0,
        max(
          0.72,
          max(
            df$TargetOnly_Stability,
            na.rm = TRUE
          ) + 0.03
        )
      ),
      ylim = c(
        0,
        1.05
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Source support vs target-domain stability",
      x = "Target-only bootstrap stability",
      y = "Weighted source PIP",
      color = NULL,
      shape = NULL
    ) +
    theme_vita(
      base_size = 11
    ) +
    ggplot2::theme(
      legend.position = if (include_legend) {
        "bottom"
      } else {
        "none"
      },
      legend.box.margin = ggplot2::margin(
        t = 4,
        b = 4
      ),
      legend.key.width = grid::unit(
        1.05,
        "cm"
      ),
      plot.margin = ggplot2::margin(
        8,
        12,
        10,
        8
      )
    )
  if (nrow(lab) > 0 &&
      has_ggrepel) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab,
        ggplot2::aes(
          x = TargetOnly_Stability,
          y = Source_PIP_Weighted,
          label = Raw_Gene
        ),
        size = 3.0,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.35,
        point.padding = 0.25,
        show.legend = FALSE
      )
  } else if (nrow(lab) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = lab,
        ggplot2::aes(
          x = TargetOnly_Stability,
          y = Source_PIP_Weighted,
          label = Raw_Gene
        ),
        vjust = -0.8,
        size = 3,
        show.legend = FALSE
      )
  }
  p
}
# ------------------------------------------------------------------------------
# Supplementary SA2 C-index dose-response
# ------------------------------------------------------------------------------
make_supp_sa2_cindex <- function() {
  if (!has_cols(
    sa2_sum,
    c(
      "Dose",
      "Avg_C_Index"
    )
  )) {
    return(NULL)
  }
  df <- sa2_sum
  df$Dose_Num <- safe_num(
    gsub(
      "%",
      "",
      df$Dose
    )
  )
  if (!("Analysis_Type" %in%
        colnames(df))) {
    df$Analysis_Type <- ifelse(
      df$Dose_Num == 100,
      "Full-source reference",
      "Repeated downsampling"
    )
  }
  df <- df %>%
    dplyr::filter(
      is.finite(Dose_Num),
      is.finite(
        safe_num(Avg_C_Index)
      )
    ) %>%
    dplyr::arrange(
      Dose_Num
    ) %>%
    dplyr::mutate(
      Type = ifelse(
        grepl(
          "Full-source",
          Analysis_Type
        ) |
          Dose_Num == 100,
        "Full-source reference",
        "Repeated downsampling"
      )
    )
  if (nrow(df) == 0) {
    return(NULL)
  }
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Dose_Num,
      y = safe_num(Avg_C_Index)
    )
  ) +
    ggplot2::geom_line(
      color = COL$orange,
      linewidth = 0.95,
      na.rm = TRUE
    ) +
    ggplot2::geom_errorbar(
      data = df %>%
        dplyr::filter(
          Type == "Repeated downsampling",
          is.finite(
            safe_num(SD_C_Index)
          )
        ),
      ggplot2::aes(
        ymin =
          safe_num(Avg_C_Index) -
          safe_num(SD_C_Index),
        ymax =
          safe_num(Avg_C_Index) +
          safe_num(SD_C_Index)
      ),
      color = COL$black,
      width = 3,
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = Type
      ),
      color = COL$orange,
      fill = "white",
      size = 3.0,
      na.rm = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Repeated downsampling" = 16,
        "Full-source reference" = 17
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(df$Dose_Num)
      )
    ) +
    ggplot2::labs(
      title = paste(
        "Source dose-response:",
        "predictive diagnostic"
      ),
      subtitle = paste(
        "Error bars are shown for repeated downsampling;",
        "100% is the full-source reference"
      ),
      x = "Source data proportion (%)",
      y = "AFT-CV C-index",
      shape = NULL
    ) +
    theme_vita(
      base_size = 12
    ) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}
# ------------------------------------------------------------------------------
# Supplementary and main-text threshold-grid plots
# ------------------------------------------------------------------------------
make_supp_threshold_grid <- function() {
  if (!has_cols(
    sa3_sum,
    c(
      "TARGET_BFDR",
      "STABILITY_CUTOFF",
      "TargetAnchored_Recall",
      "N_TargetAnchored_Selected"
    )
  )) {
    return(NULL)
  }
  df <- sa3_sum %>%
    dplyr::mutate(
      TARGET_BFDR_F = factor(
        TARGET_BFDR
      ),
      STABILITY_CUTOFF_F = factor(
        STABILITY_CUTOFF
      ),
      TargetAnchored_Recall = safe_num(
        TargetAnchored_Recall
      )
    )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = TARGET_BFDR_F,
      y = STABILITY_CUTOFF_F,
      fill = TargetAnchored_Recall
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.6
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          "n=",
          N_TargetAnchored_Selected,
          "\nR=",
          safe_round(
            TargetAnchored_Recall,
            2
          )
        )
      ),
      size = 3.5,
      color = COL$black
    ) +
    ggplot2::scale_fill_gradient(
      low = "#E8EFF6",
      high = COL$blue_dark,
      limits = c(
        0,
        1
      ),
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title = "BFDR/stability threshold grid",
      x = "BFDR tolerance",
      y = "Bootstrap stability cutoff",
      fill = "Target-anchored\nrecall"
    ) +
    theme_vita(
      base_size = 12
    )
}
make_panel_C_threshold_grid <- function() {
  p <- make_supp_threshold_grid()
  if (is.null(p)) {
    return(NULL)
  }
  p +
    ggplot2::labs(
      title = "C. BFDR/stability threshold grid",
      fill = "Target-anchored\nrecall"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(
          0.20,
          "cm"
        ),
        barheight = grid::unit(
          2.25,
          "cm"
        ),
        direction = "vertical"
      )
    ) +
    theme_vita(
      base_size = 10
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.box.spacing = grid::unit(
        1,
        "pt"
      ),
      legend.margin = ggplot2::margin(
        l = 3
      ),
      legend.title = ggplot2::element_text(
        size = 8.2,
        hjust = 0.5,
        lineheight = 0.90
      ),
      legend.text = ggplot2::element_text(
        size = 7.8
      ),
      axis.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        3,
        3,
        3,
        3
      )
    )
}
# ------------------------------------------------------------------------------
# Supplementary Cox-Snell diagnostics
# ------------------------------------------------------------------------------
make_supp_coxsnell <- function() {
  if (!has_cols(
    sa5_data,
    c(
      "Signature_Set",
      "Residual_Time",
      "Cumulative_Hazard"
    )
  )) {
    return(NULL)
  }
  df <- sa5_data %>%
    dplyr::mutate(
      Residual_Time =
        safe_num(Residual_Time),
      Cumulative_Hazard =
        safe_num(Cumulative_Hazard),
      Signature_Set =
        model_label(Signature_Set)
    ) %>%
    dplyr::filter(
      is.finite(Residual_Time),
      is.finite(Cumulative_Hazard)
    )
  if (nrow(df) == 0) {
    return(NULL)
  }
  max_axis <- max(
    c(
      df$Residual_Time,
      df$Cumulative_Hazard
    ),
    na.rm = TRUE
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Residual_Time,
      y = Cumulative_Hazard
    )
  ) +
    ggplot2::geom_step(
      color = COL$black,
      linewidth = 0.85
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = COL$dark_grey,
      linewidth = 0.75
    ) +
    ggplot2::facet_wrap(
      ~ Signature_Set,
      scales = "fixed"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(
        0,
        max_axis
      ),
      ylim = c(
        0,
        max_axis
      )
    ) +
    ggplot2::labs(
      title = "Cox-Snell residual diagnostics",
      x = "Cox-Snell residuals",
      y = "Estimated cumulative hazard"
    ) +
    theme_vita(
      base_size = 12
    )
}
# ------------------------------------------------------------------------------
# Supplementary PCA plots
# ------------------------------------------------------------------------------
make_pca_plot <- function(
    df,
    title,
    subtitle
) {
  if (!has_cols(
    df,
    c(
      "Dataset",
      "PC1",
      "PC2",
      "PC1_Var",
      "PC2_Var",
      "N_Genes"
    )
  )) {
    return(NULL)
  }
  pdat <- df %>%
    dplyr::mutate(
      PC1 = safe_num(PC1),
      PC2 = safe_num(PC2)
    ) %>%
    dplyr::filter(
      is.finite(PC1),
      is.finite(PC2)
    )
  if (nrow(pdat) == 0) {
    return(NULL)
  }
  pc1_var <- safe_first_num(
    pdat,
    "PC1_Var"
  )
  pc2_var <- safe_first_num(
    pdat,
    "PC2_Var"
  )
  n_genes <- safe_first_num(
    pdat,
    "N_Genes"
  )
  ggplot2::ggplot(
    pdat,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = Dataset
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.52,
      size = 1.8
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(
        group = Dataset
      ),
      type = "norm",
      linetype = "dashed",
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = pal_dataset
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        subtitle,
        "; candidate genes: ",
        n_genes
      ),
      x = paste0(
        "PC1 (",
        safe_round(
          pc1_var,
          2
        ),
        "% variance)"
      ),
      y = paste0(
        "PC2 (",
        safe_round(
          pc2_var,
          2
        ),
        "% variance)"
      ),
      color = NULL
    ) +
    theme_vita(
      base_size = 12
    ) +
    ggplot2::theme(
      legend.position = "top"
    )
}
# ------------------------------------------------------------------------------
# Supplementary hyperparameter grid
# ------------------------------------------------------------------------------
make_supp_hyperparameter <- function() {
  if (!has_cols(
    sa7_sum,
    c(
      "Fit_Success",
      "TargetAnchored_Recall",
      "SourceSupported_Fraction_TargetAnchored"
    )
  )) {
    return(NULL)
  }
  df <- sa7_sum %>%
    dplyr::filter(
      Fit_Success %in% TRUE
    ) %>%
    dplyr::mutate(
      TargetAnchored_Recall =
        safe_num(TargetAnchored_Recall),
      SourceSupported_Fraction_TargetAnchored =
        safe_num(
          SourceSupported_Fraction_TargetAnchored
        )
    ) %>%
    dplyr::filter(
      is.finite(TargetAnchored_Recall),
      is.finite(
        SourceSupported_Fraction_TargetAnchored
      )
    )
  if (nrow(df) == 0) {
    return(NULL)
  }
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = TargetAnchored_Recall,
      y =
        SourceSupported_Fraction_TargetAnchored
    )
  ) +
    ggplot2::geom_point(
      color = COL$black,
      alpha = 0.70,
      size = 2.0
    ) +
    ggplot2::geom_hline(
      yintercept = stats::median(
        df$SourceSupported_Fraction_TargetAnchored,
        na.rm = TRUE
      ),
      linetype = "dashed",
      color = COL$dark_grey
    ) +
    ggplot2::geom_vline(
      xintercept = stats::median(
        df$TargetAnchored_Recall,
        na.rm = TRUE
      ),
      linetype = "dashed",
      color = COL$dark_grey
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        0,
        1.05
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        1.05
      )
    ) +
    ggplot2::labs(
      title = paste(
        "Hyperparameter robustness:",
        "recall and source support"
      ),
      subtitle = paste(
        "Dashed lines show median values across",
        "successful reduced-grid cells"
      ),
      x = "Target-anchored recall",
      y = "Source-supported fraction"
    ) +
    theme_vita(
      base_size = 12
    )
}
# ------------------------------------------------------------------------------
# Main-text four-panel sensitivity figure
# ------------------------------------------------------------------------------
main_panels <- list(
  A = make_panel_A_noise(),
  B = make_panel_B_source_recovery(),
  C = make_panel_C_threshold_grid(),
  D = make_panel_D_prediction_gain()
)
main_panels <- main_panels[
  !vapply(
    main_panels,
    is.null,
    logical(1)
  )
]
if (length(main_panels) > 0) {
  panel_names <- c(
    A =
      "Figure4A_MainText_NoiseKnockin_Table",
    B =
      "Figure4B_MainText_SourceDoseResponse",
    C =
      "Figure4C_MainText_ThresholdGrid",
    D =
      "Figure4D_MainText_FixedSignaturePredictionGain"
  )
  panel_sizes <- list(
    A = c(
      5.4,
      4.1
    ),
    B = c(
      5.4,
      4.3
    ),
    C = c(
      5.4,
      4.2
    ),
    D = c(
      5.4,
      4.2
    )
  )
  for (nm in names(main_panels)) {
    size <- panel_sizes[[nm]]
    if (is.null(size)) {
      size <- c(
        5.4,
        4.2
      )
    }
    save_final_plot(
      main_panels[[nm]],
      panel_names[[nm]],
      width = size[1],
      height = size[2],
      dpi = MAIN_FIG_DPI
    )
  }
  if (has_patchwork &&
      length(main_panels) == 4) {
    composite_panels <- list(
      A = main_panels$A +
        ggplot2::theme(
          plot.margin = ggplot2::margin(
            3,
            3,
            3,
            3
          )
        ),
      B = main_panels$B +
        ggplot2::theme(
          plot.margin = ggplot2::margin(
            3,
            5,
            3,
            3
          )
        ),
      C = main_panels$C +
        ggplot2::theme(
          plot.margin = ggplot2::margin(
            3,
            3,
            3,
            3
          )
        ),
      D = main_panels$D +
        ggplot2::theme(
          plot.margin = ggplot2::margin(
            3,
            5,
            3,
            3
          )
        )
    )
    if (has_patchwork_free) {
      # Panel A has no axes or guide, so it should not inherit the space
      # reserved by Panel C's y-axis and right-side color bar or Panel B's
      # bottom legend.
      composite_panels$A <- patchwork::free(
        composite_panels$A,
        type = "panel",
        side = "lrb"
      )
    }
    main_fig <- patchwork::wrap_plots(
      composite_panels,
      ncol = 2,
      byrow = TRUE
    ) +
      patchwork::plot_layout(
        widths = c(
          1,
          1
        ),
        heights = c(
          1.02,
          0.98
        ),
        guides = "keep"
      ) +
      patchwork::plot_annotation(
        theme = ggplot2::theme(
          plot.margin = ggplot2::margin(
            2,
            2,
            2,
            2
          )
        )
      )
    save_final_plot(
      main_fig,
      "Figure4_Glioma_MainText_Sensitivity",
      width = MAIN_FIG_WIDTH_IN,
      height = MAIN_FIG_HEIGHT_IN,
      dpi = MAIN_FIG_DPI
    )
  } else if (has_patchwork) {
    main_fig <- patchwork::wrap_plots(
      main_panels,
      ncol = 2,
      byrow = TRUE
    ) +
      patchwork::plot_layout(
        guides = "keep"
      )
    save_final_plot(
      main_fig,
      "Figure4_Glioma_MainText_Sensitivity",
      width = MAIN_FIG_WIDTH_IN,
      height = MAIN_FIG_HEIGHT_IN,
      dpi = MAIN_FIG_DPI
    )
  } else {
    warning(
      "Package 'patchwork' is not installed. Main composite figure was not ",
      "assembled, but individual panels were saved."
    )
  }
} else {
  warning(
    "No main-text sensitivity panels were generated because required ",
    "inputs were unavailable."
  )
}
# ------------------------------------------------------------------------------
# Supplementary figures
# ------------------------------------------------------------------------------
supp_plots <- list(
  SA2_Cindex =
    make_supp_sa2_cindex(),
  SA3_ThresholdGrid =
    make_supp_threshold_grid(),
  SA5_CoxSnell =
    make_supp_coxsnell(),
  SA6_SourcePIP_TargetStability =
    make_panel_D_source_pip(
      include_legend = TRUE
    ),
  SA6_RawPCA =
    make_pca_plot(
      sa6_pca_raw,
      "Raw source-target feature-space separation",
      "Raw pooled candidate-gene expression space"
    ),
  SA6_WithinCohortPCA =
    make_pca_plot(
      sa6_pca_within,
      "Within-cohort standardized source-target PCA",
      paste(
        "Each cohort was gene-wise z-scored",
        "before pooled PCA"
      )
    ),
  SA7_HyperparameterGrid =
    make_supp_hyperparameter()
)
supp_plots <- supp_plots[
  !vapply(
    supp_plots,
    is.null,
    logical(1)
  )
]
if (length(supp_plots) > 0) {
  supp_sizes <- list(
    SA2_Cindex = c(
      6.6,
      4.8
    ),
    SA3_ThresholdGrid = c(
      6.2,
      5.0
    ),
    SA5_CoxSnell = c(
      10.5,
      4.8
    ),
    SA6_SourcePIP_TargetStability = c(
      8.8,
      6.4
    ),
    SA6_RawPCA = c(
      7.0,
      5.2
    ),
    SA6_WithinCohortPCA = c(
      7.0,
      5.2
    ),
    SA7_HyperparameterGrid = c(
      6.2,
      5.0
    )
  )
  for (nm in names(supp_plots)) {
    sz <- supp_sizes[[nm]]
    
    if (is.null(sz)) {
      sz <- c(
        7,
        5
      )
    }
    save_final_plot(
      supp_plots[[nm]],
      paste0(
        "Supplementary_Figure_",
        nm
      ),
      width = sz[1],
      height = sz[2],
      dpi = 600
    )
  }
}
# Optional supplementary overview figure containing all seven diagnostic panels.
# This is intended only as a compact overview for supplementary materials, not as
# the recommended main-text figure.
if (has_patchwork) {
  overview_panels <- list(
    A = main_panels$A,
    B = make_panel_B_source_recovery(),
    C = make_panel_C_threshold_grid(),
    D = make_panel_D_prediction_gain(),
    E = make_supp_coxsnell(),
    F = make_panel_D_source_pip(
      include_legend = FALSE
    ),
    G = make_supp_hyperparameter()
  )
  overview_panels <- overview_panels[
    !vapply(
      overview_panels,
      is.null,
      logical(1)
    )
  ]
  if (length(overview_panels) > 0) {
    overview_fig <- patchwork::wrap_plots(
      overview_panels,
      ncol = 2
    )
    save_final_plot(
      overview_fig,
      "Supplementary_Figure_Sensitivity_Overview",
      width = 13.5,
      height = 12.5,
      dpi = 600
    )
  }
}
# ------------------------------------------------------------------------------
# Key text summary
# ------------------------------------------------------------------------------
key_lines <- c(
  "Target-anchored glioma sensitivity analysis summary",
  "===================================================",
  "",
  "Recommended main-text Figure 4:",
  "A. Noise knock-in robustness table",
  "B. Source dose-response gene-set recovery",
  "C. BFDR/stability threshold grid",
  "D. Repeated-CV predictive gain of fixed signatures",
  "",
  "Supplementary figures generated by this script:",
  "- Source-dose C-index diagnostic (SA2)",
  "- BFDR/stability threshold grid (SA3)",
  "- Cox-Snell residual diagnostics (SA5)",
  "- Source PIP versus target-domain stability (SA6)",
  "- Raw source-target PCA (SA6)",
  "- Within-cohort standardized PCA (SA6)",
  "- Hyperparameter robustness grid (SA7)",
  "- Optional seven-panel sensitivity overview",
  "",
  "Overall summary table:",
  paste(
    capture.output(
      print(summary_table)
    ),
    collapse = "\n"
  ),
  "",
  "Interpretation notes:",
  paste(
    "- Source support is used to annotate and prioritize transferability,",
    "not to veto target-domain stability."
  ),
  paste(
    "- Target-only source-discordant genes should be interpreted as",
    "target-domain signals with weak source support."
  ),
  paste(
    "- Fixed-signature prediction remains descriptive and is not",
    "fully nested discovery-validation."
  ),
  paste(
    "- Cox-Snell residual diagnostics evaluate compatibility with",
    "the log-normal AFT assumption, not causal validity."
  ),
  paste(
    "- Source dose-response uses 30% and 60% repeated downsampling,",
    "with 100% shown as the full-source reference."
  ),
  paste(
    "- Raw source-target PCA is interpreted as a platform-separation",
    "diagnostic, not as evidence that transfer should be based on",
    "direct pooled expression overlap."
  ),
  paste(
    "- Within-cohort standardized PCA evaluates whether first-order",
    "scale differences explain the raw separation."
  ),
  paste(
    "- Source dose-response and source-target diagnostics should be",
    "interpreted as transferability diagnostics rather than",
    "confirmatory validation."
  )
)
writeLines(
  key_lines,
  final_file(
    "glioma_sensitivity_key_results.txt"
  )
)
cat(
  "\nFinal sensitivity outputs complete.\n"
)
cat(
  "Output directory:\n",
  SENS_FINAL_DIR,
  "\n",
  sep = ""
)
cat(
  "\nMain-text figure outputs:\n"
)
cat(
  "  ",
  final_file(
    "Figure4_Glioma_MainText_Sensitivity.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Figure4_Glioma_MainText_Sensitivity.pdf"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Figure4A_MainText_NoiseKnockin_Table.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Figure4B_MainText_SourceDoseResponse.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Figure4C_MainText_ThresholdGrid.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Figure4D_MainText_FixedSignaturePredictionGain.png"
  ),
  "\n",
  sep = ""
)
cat(
  "\nSupplementary figure outputs include:\n"
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA2_Cindex.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA3_ThresholdGrid.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA5_CoxSnell.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA6_SourcePIP_TargetStability.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA6_RawPCA.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA6_WithinCohortPCA.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_SA7_HyperparameterGrid.png"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "Supplementary_Figure_Sensitivity_Overview.png"
  ),
  "\n",
  sep = ""
)
cat(
  "\nKey table outputs:\n"
)
cat(
  "  ",
  final_file(
    "Table_Glioma_Sensitivity_Overall_Summary.csv"
  ),
  "\n",
  sep = ""
)
cat(
  "  ",
  final_file(
    "glioma_sensitivity_key_results.txt"
  ),
  "\n",
  sep = ""
)