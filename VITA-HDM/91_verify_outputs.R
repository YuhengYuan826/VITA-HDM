# ==============================================================================
# 91_verify_outputs.R
# Run after all analyses: Rscript 91_verify_outputs.R
# Verify that the expected publication workflow outputs exist.
#
# This script does not run or modify an analysis.
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
PROJECT_ROOT <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_ROOT)

expected_outputs <- c(
  file.path(PATHS$simulation_results, "Simulation_Raw_All.csv"),
  file.path(
    PATHS$simulation_results,
    "Final_Simulation_Summary_ForPlotting.csv"
  ),
  file.path(
    PATHS$simulation_results,
    "StructureOnly_Transfer_Sensitivity_Summary.csv"
  ),
  file.path(PATHS$case_data, "glioma_example_aligned_object.rds"),
  file.path(PATHS$discovery, "glioma_final_selection_objects.rds"),
  file.path(
    PATHS$discovery,
    "glioma_final_selected_genes_target_anchored.csv"
  ),
  file.path(
    PATHS$fixed_signature_cv,
    "glioma_fixed_signature_cindex_summary.csv"
  ),
  file.path(
    PATHS$conditional_effects,
    "glioma_overall_effects_conditional_bootstrap.csv"
  ),
  file.path(PATHS$case_summary, "glioma_example_key_results.txt"),
  file.path(
    PATHS$sensitivity_tables,
    "SA1_NoiseKnockin_MainSummary.csv"
  ),
  file.path(PATHS$sensitivity_tables, "SA2_SourceDose_Summary.csv"),
  file.path(PATHS$sensitivity_tables, "SA3_ThresholdGrid_Summary.csv"),
  file.path(
    PATHS$sensitivity_tables,
    "SA4_FixedSignature_Prediction_Summary.csv"
  ),
  file.path(
    PATHS$sensitivity_tables,
    "SA5_CoxSnell_Residual_Summary.csv"
  ),
  file.path(
    PATHS$sensitivity_tables,
    "SA6_SourcePIP_TargetStability_Correlation.csv"
  ),
  file.path(
    PATHS$sensitivity_tables,
    "SA7_HyperparameterGrid_Summary.csv"
  ),
  file.path(
    PATHS$sensitivity_tables,
    "SA8_TargetAnchoring_RuleSummary.csv"
  ),
  file.path(
    PATHS$sensitivity_final,
    "Figure4_Glioma_MainText_Sensitivity.png"
  ),
  file.path(
    PATHS$sensitivity_final,
    "Table_Glioma_Sensitivity_Overall_Summary.csv"
  ),
  file.path(
    PATHS$sensitivity_final,
    "glioma_sensitivity_key_results.txt"
  )
)

status <- data.frame(
  Path = expected_outputs,
  Exists = file.exists(expected_outputs),
  stringsAsFactors = FALSE
)
print(status, row.names = FALSE)

missing_outputs <- status$Path[!status$Exists]
if (length(missing_outputs) > 0) {
  stop(
    "Expected outputs are missing:\n",
    paste(missing_outputs, collapse = "\n"),
    call. = FALSE
  )
}

private_pca_files <- c(
  "SA6_PCA_SourceTarget_Data.csv",
  "SA6_PCA_SourceTarget_RawPooled_Data.csv",
  "SA6_PCA_SourceTarget_WithinCohortZ_Data.csv"
)
misplaced_private <- file.exists(file.path(
  PATHS$sensitivity_final,
  private_pca_files
))

if (any(misplaced_private)) {
  stop(
    "Participant-level SA6 PCA files were found in final_outputs: ",
    paste(private_pca_files[misplaced_private], collapse = ", "),
    call. = FALSE
  )
}

message("Expected outputs are complete and private SA6 tables are isolated.")
