# ==============================================================================
# 12_summarise_glioma_results.R
# Run after files 10 and 11: Rscript 12_summarise_glioma_results.R
# Summarise glioma case-study outputs for the manuscript and supplement.
# Analysis specification:
#   - Summarises data preparation, discovery, CV, and effect outputs.
#   - Uses the target-anchored reported set as the main real-data gene set.
#   - Separately reports:
#       1. target-anchored reported set
#       2. transfer-supported subset
#       3. target-only source-discordant genes
#       4. target-only stable genes
#       5. fixed-signature predictive assessment
#       6. conditional bootstrap effects
# Inputs:
#   results/case_study/
#     data/
#     discovery/
#     fixed_signature_cv/
#     conditional_effects/
# Outputs:
#   results/case_study/summary/
#     Table_RealData_CohortSummary.csv
#     Table_RealData_FinalReportedSet.csv
#     Table_RealData_FinalSignature.csv
#     Table_RealData_TransferSupportedSet.csv
#     Table_RealData_TargetOnlySourceDiscordant.csv
#     Table_RealData_TargetOnlyStableGenes.csv
#     Table_RealData_PredictiveAssessment.csv
#     Table_RealData_PredictiveContrasts.csv
#     Table_RealData_ConditionalEffects.csv
#     Table_Supplement_GeneLevelConditionalEffects.csv
#     glioma_example_key_results.txt
# No sink() logging is used.
# ==============================================================================
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
PROJECT_DIR <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_DIR)
suppressPackageStartupMessages({
  library(dplyr)
})
RESULTS_ROOT <- PATHS$case_results
DATA_DIR <- PATHS$case_data
DISC_DIR <- PATHS$discovery
CV_DIR <- PATHS$fixed_signature_cv
EFF_DIR <- PATHS$conditional_effects
OUT_DIR <- ensure_directory(PATHS$case_summary)
# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
first_existing_file <- function(paths, label = "file", required = TRUE) {
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Cannot find required ", label, ". Tried:\n",
      paste(paths, collapse = "\n")
    )
  }
  NA_character_
}
read_csv_if_exists <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(data.frame())
  }
  read.csv(path, stringsAsFactors = FALSE)
}
round_numeric_cols <- function(df, digits = 4) {
  if (is.null(df) || nrow(df) == 0) return(df)
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (cc in num_cols) df[[cc]] <- round(df[[cc]], digits)
  df
}
safe_pick_cols <- function(df, cols) {
  cols <- intersect(cols, colnames(df))
  df[, cols, drop = FALSE]
}
safe_count_true <- function(df, col) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% colnames(df))) {
    return(NA_integer_)
  }
  sum(df[[col]] %in% TRUE, na.rm = TRUE)
}
safe_fraction_true <- function(df, col) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% colnames(df))) {
    return(NA_real_)
  }
  mean(df[[col]] %in% TRUE, na.rm = TRUE)
}
safe_value <- function(df, row_expr, col) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% colnames(df))) {
    return(NA)
  }
  idx <- which(row_expr)
  if (length(idx) == 0) return(NA)
  df[[col]][idx[1]]
}
# ------------------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------------------
cohort_summary_file <- file.path(DATA_DIR, "glioma_example_cohort_summary.csv")
gene_overlap_file <- file.path(DATA_DIR, "glioma_example_gene_overlap_summary.csv")
source_ratio_file <- file.path(DATA_DIR, "glioma_example_source_target_ratio.csv")
expression_summary_file <- file.path(DATA_DIR, "glioma_example_expression_summary.csv")
tcga_idh_diag_file <- file.path(DATA_DIR, "glioma_example_tcga_idh_diagnostics.csv")
target_anchored_file <- first_existing_file(
  c(
    file.path(DISC_DIR, "glioma_final_selected_genes_target_anchored.csv"),
    file.path(DISC_DIR, "glioma_final_selected_genes_manuscript.csv")
  ),
  label = "target-anchored reported-set file",
  required = TRUE
)
transfer_supported_file <- first_existing_file(
  c(file.path(DISC_DIR, "glioma_final_selected_genes_transfer_supported.csv")),
  label = "transfer-supported set file",
  required = FALSE
)
target_discordant_file <- first_existing_file(
  c(file.path(DISC_DIR, "glioma_targetonly_source_discordant_genes.csv")),
  label = "target-only source-discordant gene file",
  required = FALSE
)
targetonly_file <- first_existing_file(
  c(file.path(DISC_DIR, "glioma_targetonly_stable_genes.csv")),
  label = "target-only stable gene file",
  required = FALSE
)
selection_summary_file <- file.path(DISC_DIR, "glioma_selection_summary.csv")
selection_report_file <- file.path(DISC_DIR, "glioma_bootstrap_selection_report.csv")
source_weights_file <- file.path(DISC_DIR, "glioma_source_weights_full.csv")
rule_file <- file.path(DISC_DIR, "glioma_manuscript_signature_rule.csv")
cindex_summary_file <- file.path(CV_DIR, "glioma_fixed_signature_cindex_summary.csv")
cindex_contrasts_file <- first_existing_file(
  c(file.path(CV_DIR, "glioma_fixed_signature_cindex_pairwise_contrasts.csv")),
  label = "fixed-signature pairwise contrast file",
  required = FALSE
)
overall_effects_file <- file.path(EFF_DIR, "glioma_overall_effects_conditional_bootstrap.csv")
gene_effects_file <- file.path(EFF_DIR, "glioma_gene_mediation_effects_conditional_bootstrap.csv")
effect_meta_file <- file.path(EFF_DIR, "glioma_conditional_effect_run_metadata.csv")
required_files <- c(
  cohort_summary_file,
  gene_overlap_file,
  source_ratio_file,
  target_anchored_file,
  selection_summary_file,
  source_weights_file,
  cindex_summary_file,
  overall_effects_file,
  gene_effects_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required files. Run data preparation, discovery, fixed-signature ",
    "CV, and conditional-effect scripts first:\n",
    paste(missing_files, collapse = "\n")
  )
}
# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
cohort_summary <- read.csv(cohort_summary_file, stringsAsFactors = FALSE)
gene_overlap <- read.csv(gene_overlap_file, stringsAsFactors = FALSE)
source_ratio <- read.csv(source_ratio_file, stringsAsFactors = FALSE)
expression_summary <- read_csv_if_exists(expression_summary_file)
tcga_idh_diag <- read_csv_if_exists(tcga_idh_diag_file)
selected_target_anchored <- read.csv(target_anchored_file, stringsAsFactors = FALSE)
selected_transfer_supported <- read_csv_if_exists(transfer_supported_file)
selected_target_discordant <- read_csv_if_exists(target_discordant_file)
targetonly_stable <- read_csv_if_exists(targetonly_file)
selection_summary <- read.csv(selection_summary_file, stringsAsFactors = FALSE)
selection_report <- read_csv_if_exists(selection_report_file)
source_weights <- read.csv(source_weights_file, stringsAsFactors = FALSE)
rule_df <- read_csv_if_exists(rule_file)
cindex_summary <- read.csv(cindex_summary_file, stringsAsFactors = FALSE)
cindex_contrasts <- read_csv_if_exists(cindex_contrasts_file)
overall_effects <- read.csv(overall_effects_file, stringsAsFactors = FALSE)
gene_effects <- read.csv(gene_effects_file, stringsAsFactors = FALSE)
effect_meta <- read_csv_if_exists(effect_meta_file)
# ------------------------------------------------------------------------------
# Ensure Raw_Gene columns
# ------------------------------------------------------------------------------
ensure_raw_gene <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!("Raw_Gene" %in% colnames(df))) {
    if ("Gene" %in% colnames(df)) {
      df$Raw_Gene <- sub("^M_", "", df$Gene)
    }
  }
  df
}
selected_target_anchored <- ensure_raw_gene(selected_target_anchored)
selected_transfer_supported <- ensure_raw_gene(selected_transfer_supported)
selected_target_discordant <- ensure_raw_gene(selected_target_discordant)
targetonly_stable <- ensure_raw_gene(targetonly_stable)
gene_effects <- ensure_raw_gene(gene_effects)
# ------------------------------------------------------------------------------
# Manuscript-ready cohort summary
# ------------------------------------------------------------------------------
table_cohort <- cohort_summary %>%
  transmute(
    Cohort,
    N = N_Complete_Basic,
    Events,
    Censored,
    Event_Rate = round(Event_Rate, 3),
    IDH_Mutant,
    IDH_Wildtype,
    N_With_Grade,
    Median_OS_Time = round(Median_OS_Time, 2),
    N_Expression_Genes
  )
write.csv(
  table_cohort,
  file.path(OUT_DIR, "Table_RealData_CohortSummary.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Expression and IDH diagnostic tables, optional
# ------------------------------------------------------------------------------
if (nrow(expression_summary) > 0) {
  table_expression <- expression_summary %>%
    mutate(
      Min = round(Min, 4),
      Q25 = round(Q25, 4),
      Median = round(Median, 4),
      Q75 = round(Q75, 4),
      Q99 = round(Q99, 4),
      Max = round(Max, 4)
    )
  write.csv(
    table_expression,
    file.path(OUT_DIR, "Table_Supplement_ExpressionSummary.csv"),
    row.names = FALSE
  )
}
if (nrow(tcga_idh_diag) > 0) {
  write.csv(
    tcga_idh_diag,
    file.path(OUT_DIR, "Table_Supplement_TCGA_IDH_Diagnostics.csv"),
    row.names = FALSE
  )
}
# ------------------------------------------------------------------------------
# Final target-anchored reported set table
# ------------------------------------------------------------------------------
signature_cols <- c(
  "Raw_Gene",
  "Category",
  "Final_Target_PIP_VITA",
  "Final_Target_Beta_VITA",
  "Final_Target_PIP_TargetOnly",
  "Final_Target_Beta_TargetOnly",
  "Source_PIP_Weighted",
  "Source_Beta_Weighted",
  "VITA_Stability",
  "TargetOnly_Stability",
  "Source_Supported",
  "VITA_Extended",
  "VITA_Core",
  "TargetOnly_Stable",
  "SourceSupported_Stable",
  "Target_Enhanced",
  "TargetOnly_SourceDiscordant",
  "TransferSupported_Signature",
  "TargetAnchored_Reported",
  "FullData_Selected_Source",
  "FullData_Selected_TargetOnly",
  "FullData_Selected_VITA"
)
table_target_anchored <- selected_target_anchored %>%
  safe_pick_cols(signature_cols) %>%
  round_numeric_cols(digits = 4)
if ("VITA_Stability" %in% colnames(table_target_anchored)) {
  table_target_anchored <- table_target_anchored %>%
    arrange(desc(VITA_Stability), desc(TargetOnly_Stability), Raw_Gene)
}
write.csv(
  table_target_anchored,
  file.path(OUT_DIR, "Table_RealData_FinalReportedSet.csv"),
  row.names = FALSE
)
# Backward-compatible filename.
write.csv(
  table_target_anchored,
  file.path(OUT_DIR, "Table_RealData_FinalSignature.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Transfer-supported subset table
# ------------------------------------------------------------------------------
table_transfer_supported <- selected_transfer_supported %>%
  safe_pick_cols(signature_cols) %>%
  round_numeric_cols(digits = 4)
if (nrow(table_transfer_supported) > 0 &&
    "VITA_Stability" %in% colnames(table_transfer_supported)) {
  table_transfer_supported <- table_transfer_supported %>%
    arrange(desc(VITA_Stability), Raw_Gene)
}
write.csv(
  table_transfer_supported,
  file.path(OUT_DIR, "Table_RealData_TransferSupportedSet.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Target-only source-discordant genes table
# ------------------------------------------------------------------------------
table_target_discordant <- selected_target_discordant %>%
  safe_pick_cols(signature_cols) %>%
  round_numeric_cols(digits = 4)
if (nrow(table_target_discordant) > 0 &&
    "TargetOnly_Stability" %in% colnames(table_target_discordant)) {
  table_target_discordant <- table_target_discordant %>%
    arrange(desc(TargetOnly_Stability), Raw_Gene)
}
write.csv(
  table_target_discordant,
  file.path(OUT_DIR, "Table_RealData_TargetOnlySourceDiscordant.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Target-only stable genes table
# ------------------------------------------------------------------------------
table_targetonly <- targetonly_stable %>%
  safe_pick_cols(signature_cols) %>%
  round_numeric_cols(digits = 4)
if (nrow(table_targetonly) > 0 &&
    "TargetOnly_Stability" %in% colnames(table_targetonly)) {
  table_targetonly <- table_targetonly %>%
    arrange(desc(TargetOnly_Stability), Raw_Gene)
}
write.csv(
  table_targetonly,
  file.path(OUT_DIR, "Table_RealData_TargetOnlyStableGenes.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Predictive assessment table
# ------------------------------------------------------------------------------
table_pred <- cindex_summary %>%
  transmute(
    Model,
    Training_Domain,
    N_Genes,
    RepeatedCV_N,
    Mean_C_Index = round(RepeatedCV_Mean_C_Index, 4),
    SD_C_Index = round(RepeatedCV_SD_C_Index, 4),
    Mean_Delta_vs_Clinical = round(Mean_Delta_vs_Clinical, 4),
    Positive_Delta_N,
    Apparent_C_Index = round(Apparent_C_Index, 4)
  ) %>%
  arrange(desc(Mean_C_Index))
write.csv(
  table_pred,
  file.path(OUT_DIR, "Table_RealData_PredictiveAssessment.csv"),
  row.names = FALSE
)
if (nrow(cindex_contrasts) > 0) {
  table_pred_contrasts <- cindex_contrasts %>%
    mutate(
      Mean_Difference = round(Mean_Difference, 4),
      SD_Difference = round(SD_Difference, 4)
    )
  write.csv(
    table_pred_contrasts,
    file.path(OUT_DIR, "Table_RealData_PredictiveContrasts.csv"),
    row.names = FALSE
  )
} else {
  table_pred_contrasts <- data.frame()
  write.csv(
    table_pred_contrasts,
    file.path(OUT_DIR, "Table_RealData_PredictiveContrasts.csv"),
    row.names = FALSE
  )
}
# ------------------------------------------------------------------------------
# Conditional effect summary table
# ------------------------------------------------------------------------------
table_effects <- overall_effects %>%
  transmute(
    Effect_Set,
    N_Genes,
    Parameter,
    Estimate = round(Estimate, 4),
    Bootstrap_Mean = round(Bootstrap_Mean, 4),
    Bootstrap_SD = round(Bootstrap_SD, 4),
    CI_95_Low = round(CI_95_Low, 4),
    CI_95_Upp = round(CI_95_Upp, 4),
    CI_Crosses_Zero,
    N_Boot_Used
  )
write.csv(
  table_effects,
  file.path(OUT_DIR, "Table_RealData_ConditionalEffects.csv"),
  row.names = FALSE
)
gene_effect_cols <- c(
  "Raw_Gene",
  "Category",
  "In_TargetAnchored_Reported_Set",
  "In_TransferSupported_Set",
  "In_Source_Supported_Subset",
  "In_Target_Enhanced_Subset",
  "In_TargetOnly_Stable_Signature",
  "In_TargetOnly_SourceDiscordant_Subset",
  "Alpha_D_to_M",
  "Beta_M_to_Y",
  "PIP_Full",
  "IDE_Point",
  "IDE_Boot_Mean",
  "IDE_Boot_SD",
  "IDE_95_Low",
  "IDE_95_Upp",
  "IDE_Crosses_Zero",
  "IDE_Sign_Stability",
  "N_Boot_Used"
)
table_gene_effects <- gene_effects %>%
  safe_pick_cols(gene_effect_cols) %>%
  round_numeric_cols(digits = 4)
write.csv(
  table_gene_effects,
  file.path(OUT_DIR, "Table_Supplement_GeneLevelConditionalEffects.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Key-result extraction
# ------------------------------------------------------------------------------
n_reported <- nrow(selected_target_anchored)
n_transfer <- nrow(selected_transfer_supported)
n_target_discordant <- nrow(selected_target_discordant)
n_targetonly <- nrow(targetonly_stable)
source_supported_n <- safe_count_true(selected_target_anchored, "Source_Supported")
source_supported_frac <- safe_fraction_true(selected_target_anchored, "Source_Supported")
core_n <- safe_count_true(selected_target_anchored, "VITA_Core")
transfer_source_supported_n <- safe_count_true(
  selected_transfer_supported,
  "Source_Supported"
)
target_anchored_genes_text <- if (nrow(table_target_anchored) > 0 &&
                                  "Raw_Gene" %in% colnames(table_target_anchored)) {
  paste(table_target_anchored$Raw_Gene, collapse = ", ")
} else {
  "None"
}
transfer_genes_text <- if (nrow(table_transfer_supported) > 0 &&
                           "Raw_Gene" %in% colnames(table_transfer_supported)) {
  paste(table_transfer_supported$Raw_Gene, collapse = ", ")
} else {
  "None"
}
discordant_genes_text <- if (nrow(table_target_discordant) > 0 &&
                             "Raw_Gene" %in% colnames(table_target_discordant)) {
  paste(table_target_discordant$Raw_Gene, collapse = ", ")
} else {
  "None"
}
targetonly_genes_text <- if (nrow(table_targetonly) > 0 &&
                             "Raw_Gene" %in% colnames(table_targetonly)) {
  paste(table_targetonly$Raw_Gene, collapse = ", ")
} else {
  "None"
}
clin_row <- table_pred[table_pred$Model == "Clinical-Only", , drop = FALSE]
anchored_row <- table_pred[
  grepl("TargetAnchored_Reported_Set", table_pred$Model),
  ,
  drop = FALSE
]
transfer_row <- table_pred[
  grepl("TransferSupported_Set", table_pred$Model),
  ,
  drop = FALSE
]
source_row <- table_pred[
  grepl("^Source-Only", table_pred$Model),
  ,
  drop = FALSE
]
targetonly_row <- table_pred[
  grepl("TargetOnly_Stable", table_pred$Model),
  ,
  drop = FALSE
]
main_effects <- table_effects[
  table_effects$Effect_Set == "TargetAnchored_Reported_Set",
  ,
  drop = FALSE
]
transfer_effects <- table_effects[
  table_effects$Effect_Set == "TransferSupported_Set",
  ,
  drop = FALSE
]
discordant_effects <- table_effects[
  table_effects$Effect_Set == "TargetOnly_SourceDiscordant_Subset",
  ,
  drop = FALSE
]
ide_row <- main_effects[
  main_effects$Parameter == "Total_Indirect_Effect",
  ,
  drop = FALSE
]
de_row <- main_effects[
  main_effects$Parameter == "Direct_Effect",
  ,
  drop = FALSE
]
te_row <- main_effects[
  main_effects$Parameter == "Total_Effect",
  ,
  drop = FALSE
]
mp_row <- main_effects[
  main_effects$Parameter == "Mediation_Proportion",
  ,
  drop = FALSE
]
# ------------------------------------------------------------------------------
# Key-result text file for manual manuscript drafting
# ------------------------------------------------------------------------------
key_lines <- c(
  "Final cleaned glioma VITA-HDM example key results",
  "==================================================",
  "",
  "Analysis rule:",
  "The final reported set is target-anchored. Target-domain stable signals are retained and annotated. Source evidence is used to prioritize transferability but is not used as a veto against target-only stable genes.",
  "",
  "Gene universe:",
  paste0("Definition: ", gene_overlap$Gene_Universe_Definition[1]),
  paste0("Main common genes: ", gene_overlap$Main_Common_Genes_Target_Source[1]),
  "",
  "Cohort summary:",
  paste(capture.output(print(table_cohort)), collapse = "\n"),
  "",
  "Source-target ratio:",
  paste(capture.output(print(source_ratio)), collapse = "\n"),
  "",
  "Target-anchored final reported set:",
  paste0("Size: ", n_reported),
  paste0("Genes: ", target_anchored_genes_text),
  paste0("Source-supported genes: ", source_supported_n,
         " (fraction = ", round(source_supported_frac, 3), ")"),
  paste0("High-stability transfer core genes: ", core_n),
  "",
  "Transfer-supported subset:",
  paste0("Size: ", n_transfer),
  paste0("Genes: ", transfer_genes_text),
  paste0("Source-supported genes in transfer-supported subset: ",
         transfer_source_supported_n),
  "",
  "Target-only source-discordant genes:",
  paste0("Size: ", n_target_discordant),
  paste0("Genes: ", discordant_genes_text),
  "",
  "Target-only stable genes:",
  paste0("Size: ", n_targetonly),
  paste0("Genes: ", targetonly_genes_text),
  "",
  "Predictive assessment:",
  paste(capture.output(print(table_pred)), collapse = "\n"),
  "",
  "Predictive pairwise contrasts:",
  if (nrow(table_pred_contrasts) > 0) {
    paste(capture.output(print(table_pred_contrasts)), collapse = "\n")
  } else {
    "No pairwise contrast table was available."
  },
  "",
  "Conditional bootstrap effect summary for target-anchored reported set:",
  if (nrow(main_effects) > 0) {
    paste(capture.output(print(main_effects)), collapse = "\n")
  } else {
    "No TargetAnchored_Reported_Set rows found in conditional effects."
  },
  "",
  "Conditional bootstrap effect summary for transfer-supported set:",
  if (nrow(transfer_effects) > 0) {
    paste(capture.output(print(transfer_effects)), collapse = "\n")
  } else {
    "No TransferSupported_Set rows found in conditional effects."
  },
  "",
  "Conditional bootstrap effect summary for target-only source-discordant subset:",
  if (nrow(discordant_effects) > 0) {
    paste(capture.output(print(discordant_effects)), collapse = "\n")
  } else {
    "No TargetOnly_SourceDiscordant_Subset rows found in conditional effects."
  },
  "",
  "Interpretation notes:",
  "- Fixed-signature CV is descriptive and not fully nested discovery-validation.",
  "- Conditional bootstrap is for fixed selected mediator sets and is not full post-selection inference.",
  "- In this single-source main example, LPML source weight degenerates to 1; LPML is retained as a transferability diagnostic rather than as a multi-source weighting mechanism.",
  "- The target-anchored rule prevents source-domain evidence from acting as a veto against target-only stable genes.",
  "- Genes in the target-only source-discordant category should be interpreted as target-domain signals with weak source support, not as biologically irrelevant genes."
)
writeLines(
  key_lines,
  file.path(OUT_DIR, "glioma_example_key_results.txt")
)
# ------------------------------------------------------------------------------
# Compact machine-readable metadata
# ------------------------------------------------------------------------------
summary_meta <- data.frame(
  Item = c(
    "N_TargetAnchored_Reported",
    "N_TransferSupported",
    "N_TargetOnly_SourceDiscordant",
    "N_TargetOnly_Stable",
    "N_SourceSupported_In_TargetAnchored",
    "Frac_SourceSupported_In_TargetAnchored",
    "N_HighStability_Core_In_TargetAnchored",
    "TargetAnchored_Genes",
    "TransferSupported_Genes",
    "TargetOnly_SourceDiscordant_Genes",
    "TargetOnly_Stable_Genes"
  ),
  Value = c(
    n_reported,
    n_transfer,
    n_target_discordant,
    n_targetonly,
    source_supported_n,
    round(source_supported_frac, 4),
    core_n,
    target_anchored_genes_text,
    transfer_genes_text,
    discordant_genes_text,
    targetonly_genes_text
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_meta,
  file.path(OUT_DIR, "Table_RealData_TargetAnchoredSummaryMetadata.csv"),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Console output
# ------------------------------------------------------------------------------
cat("\nSummary complete.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")
cat("\nKey outputs:\n")
cat("  ", file.path(OUT_DIR, "Table_RealData_CohortSummary.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_FinalReportedSet.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_TransferSupportedSet.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_TargetOnlySourceDiscordant.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_TargetOnlyStableGenes.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_PredictiveAssessment.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_PredictiveContrasts.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_RealData_ConditionalEffects.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "Table_Supplement_GeneLevelConditionalEffects.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "glioma_example_key_results.txt"), "\n", sep = "")
cat("\nTarget-anchored reported genes:\n")
cat("  ", target_anchored_genes_text, "\n", sep = "")
cat("\nTarget-only source-discordant genes:\n")
cat("  ", discordant_genes_text, "\n", sep = "")