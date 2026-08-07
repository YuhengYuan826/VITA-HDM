# ==============================================================================
# 07_prepare_glioma_data.R
# Run: Rscript 07_prepare_glioma_data.R
# Purpose:
#   Prepare the two-cohort glioma case study for VITA-HDM.
# Analysis design:
#   Target:  CGGA_mRNA_array_301
#   Source:  TCGA_LGG
#   Exposure: IDH mutation status
#   Outcome: overall survival
#   Covariates: Age + Gender
# Required input files under the repository data directory:
#   data/CGGA/
#     CGGA.mRNA_array_301_clinical.20200506.txt
#     CGGA.mRNA_array_301_gene_level.20200506.txt
#   data/TCGA/
#     TCGA-LGG.clinical.tsv
#     TCGA-LGG.somaticmutation_wxs.tsv
#     TCGA-LGG.star_fpkm.tsv
#     TCGA-LGG.survival.tsv
# Outputs are written only under:
#   results/case_study/data/
# Important:
#   - Only the target and source used in the paper define the gene universe.
#   - No CGGA325, CGGA693, or TCGA_GBM information is used here.
#   - This avoids old-result contamination and avoids restricting the final gene
#     universe by cohorts not used in the main example.
#   - No sink() logging is used; progress and errors are printed directly in R.
#   - Expression and IDH handling:
#     * expression matrices are not log2-transformed here;
#     * TCGA IDH status is inferred primarily from somatic mutation calls;
#     * clinical IDH, if present, is retained only as a diagnostic column.
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
  library(data.table)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
RESULTS_ROOT <- PATHS$case_results
OUT_DIR <- ensure_directory(PATHS$case_data)
CGGA_DIR <- PATHS$cgga_data
TCGA_DIR <- PATHS$tcga_data
CGGA_ARRAY_CLIN <- file.path(CGGA_DIR, "CGGA.mRNA_array_301_clinical.20200506.txt")
CGGA_ARRAY_EXPR <- file.path(CGGA_DIR, "CGGA.mRNA_array_301_gene_level.20200506.txt")
TCGA_LGG_CLIN <- file.path(TCGA_DIR, "TCGA-LGG.clinical.tsv")
TCGA_LGG_SURV <- file.path(TCGA_DIR, "TCGA-LGG.survival.tsv")
TCGA_LGG_MUT  <- file.path(TCGA_DIR, "TCGA-LGG.somaticmutation_wxs.tsv")
TCGA_LGG_EXPR <- file.path(TCGA_DIR, "TCGA-LGG.star_fpkm.tsv")
TARGET_COHORT <- "CGGA_mRNA_array_301"
SOURCE_COHORT <- "TCGA_LGG"
# ------------------------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------------------------
read_any_table <- function(path) {
  if (!file.exists(path)) stop("Cannot find file: ", path)
  data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
}
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clean_id <- function(x) stringr::str_trim(as.character(x))
standardize_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  toupper(x)
}
normalize_tcga_barcode <- function(x) {
  gsub("\\.", "-", as.character(x))
}
tcga_patient_id <- function(x) substr(normalize_tcga_barcode(x), 1, 12)
tcga_sample_short <- function(x) substr(normalize_tcga_barcode(x), 1, 16)
parse_gender <- function(x) {
  x0 <- tolower(stringr::str_trim(as.character(x)))
  out <- rep(NA_real_, length(x0))
  out[x0 %in% c("male", "m", "1")] <- 1
  out[x0 %in% c("female", "f", "0")] <- 0
  out
}
parse_idh <- function(x) {
  x0 <- tolower(stringr::str_trim(as.character(x)))
  out <- rep(NA_real_, length(x0))
  out[grepl("mut", x0)] <- 1
  out[grepl("wild|wt", x0)] <- 0
  out[x0 %in% c("1", "idh-mut", "idh mutant", "mutant", "mut")] <- 1
  out[x0 %in% c("0", "idh-wt", "idh wildtype", "wildtype", "wt")] <- 0
  out
}
parse_grade <- function(x) {
  x0 <- toupper(stringr::str_trim(as.character(x)))
  out <- rep(NA_real_, length(x0))
  out[grepl("II", x0) | grepl("\\b2\\b", x0)] <- 2
  out[grepl("III", x0) | grepl("\\b3\\b", x0)] <- 3
  out[grepl("IV", x0) | grepl("\\b4\\b", x0)] <- 4
  out
}
first_existing_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) stop("Cannot find ", label, ". Candidate names: ", paste(candidates, collapse = ", "))
    return(NULL)
  }
  hit[1]
}
# ------------------------------------------------------------------------------
# Expression helpers
# ------------------------------------------------------------------------------
collapse_duplicate_genes <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  keep <- !is.na(rownames(expr_mat)) &
    rownames(expr_mat) != "" &
    rownames(expr_mat) != "NA"
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (any(duplicated(rownames(expr_mat)))) {
    expr_df <- data.frame(Gene = rownames(expr_mat), expr_mat, check.names = FALSE)
    expr_df <- expr_df %>%
      group_by(Gene) %>%
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
    rownames(expr_mat) <- expr_df$Gene
  }
  expr_mat
}
auto_log2_transform <- function(expr_mat, dataset_name) {
  # No-log2 analysis specification:
  # Keep expression values exactly as read from the input files after gene-symbol
  # mapping and duplicate-gene collapsing. Downstream train-only standardization is
  # still performed inside model-design construction.
  expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  message(dataset_name, ": no log2(x + 1) transformation was applied.")
  attr(expr_mat, "log2_transformed_by_script") <- FALSE
  attr(expr_mat, "expression_transform_mode") <- "none_no_log2_audit"
  expr_mat
}
summarise_expr <- function(expr_mat, dataset_name) {
  vals <- as.numeric(expr_mat)
  vals <- vals[is.finite(vals)]
  data.frame(
    Cohort = dataset_name,
    N_Genes = nrow(expr_mat),
    N_Samples = ncol(expr_mat),
    Min = min(vals, na.rm = TRUE),
    Q25 = as.numeric(quantile(vals, 0.25, na.rm = TRUE)),
    Median = median(vals, na.rm = TRUE),
    Q75 = as.numeric(quantile(vals, 0.75, na.rm = TRUE)),
    Q99 = as.numeric(quantile(vals, 0.99, na.rm = TRUE)),
    Max = max(vals, na.rm = TRUE),
    Log2_Transformed_By_Script = isTRUE(attr(expr_mat, "log2_transformed_by_script")),
    stringsAsFactors = FALSE
  )
}
parse_cgga_expr <- function(path, dataset_name) {
  expr_df <- read_any_table(path)
  cn <- colnames(expr_df)
  gene_col_candidates <- c(
    "Gene_Name", "Gene", "gene", "Symbol", "symbol",
    "GeneSymbol", "gene_symbol", "gene_name", "V1", "ID"
  )
  gene_col <- intersect(gene_col_candidates, cn)
  gene_col <- if (length(gene_col) > 0) gene_col[1] else cn[1]
  genes <- standardize_gene_symbol(expr_df[[gene_col]])
  sample_cols <- setdiff(cn, gene_col)
  numeric_rate <- sapply(sample_cols, function(cc) {
    vals <- suppressWarnings(as.numeric(expr_df[[cc]]))
    mean(is.finite(vals))
  })
  sample_cols <- sample_cols[numeric_rate > 0.5]
  if (length(sample_cols) == 0) stop(dataset_name, ": no numeric expression sample columns detected.")
  expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- genes
  colnames(expr_mat) <- clean_id(colnames(expr_mat))
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_mat <- auto_log2_transform(expr_mat, dataset_name)
  message("\nCGGA expression parsed: ", dataset_name)
  message("  Genes: ", nrow(expr_mat))
  message("  Samples: ", ncol(expr_mat))
  expr_mat
}
parse_tcga_expr_ensembl_to_symbol <- function(path, dataset_name) {
  expr_df <- read_any_table(path)
  gene_col <- first_existing_col(
    expr_df,
    candidates = c("Ensembl_ID", "ensembl_gene_id", "gene_id", "Gene", "gene", colnames(expr_df)[1]),
    required = TRUE,
    label = paste0(dataset_name, " Ensembl/gene ID column")
  )
  ensembl_raw <- as.character(expr_df[[gene_col]])
  ensembl_clean <- sub("\\..*$", "", ensembl_raw)
  gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  gene_symbol <- standardize_gene_symbol(gene_symbol)
  sample_cols <- setdiff(colnames(expr_df), gene_col)
  expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- gene_symbol
  keep <- !is.na(rownames(expr_mat)) &
    rownames(expr_mat) != "" &
    rownames(expr_mat) != "NA"
  expr_mat <- expr_mat[keep, , drop = FALSE]
  colnames(expr_mat) <- normalize_tcga_barcode(colnames(expr_mat))
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_mat <- auto_log2_transform(expr_mat, dataset_name)
  message("\nTCGA expression parsed: ", dataset_name)
  message("  Raw rows: ", length(ensembl_raw))
  message("  Mapped unique symbols: ", nrow(expr_mat))
  message("  Samples: ", ncol(expr_mat))
  expr_mat
}
choose_primary_tcga_samples_keep_sample_id <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  sample_ids <- normalize_tcga_barcode(colnames(expr_mat))
  df <- data.frame(
    sample_raw = colnames(expr_mat),
    Sample_ID = tcga_sample_short(sample_ids),
    Patient_ID = tcga_patient_id(sample_ids),
    sample_type = substr(sample_ids, 14, 15),
    stringsAsFactors = FALSE
  )
  df$priority <- ifelse(df$sample_type == "01", 1L, 2L)
  df <- df[order(df$Patient_ID, df$priority, df$Sample_ID), , drop = FALSE]
  df <- df[!duplicated(df$Patient_ID), , drop = FALSE]
  if (nrow(df) == 0) stop("No TCGA samples remain after primary-sample selection.")
  expr_sub <- expr_mat[, df$sample_raw, drop = FALSE]
  colnames(expr_sub) <- df$Sample_ID
  if (any(duplicated(colnames(expr_sub)))) stop("Duplicated TCGA sample IDs remain.")
  expr_sub
}
# ------------------------------------------------------------------------------
# Clinical helpers
# ------------------------------------------------------------------------------
make_cgga_clin <- function(df, cohort_name) {
  data.frame(
    Sample_ID = clean_id(df[["CGGA_ID"]]),
    Patient_ID = clean_id(df[["CGGA_ID"]]),
    Cohort = cohort_name,
    OS_time = as_num(df[["OS"]]),
    OS_status = as_num(df[["Censor (alive=0; dead=1)"]]),
    IDH = parse_idh(df[["IDH_mutation_status"]]),
    Age = as_num(df[["Age"]]),
    Gender = parse_gender(df[["Gender"]]),
    Grade = if ("Grade" %in% colnames(df)) parse_grade(df[["Grade"]]) else NA_real_,
    Histology = if ("Histology" %in% colnames(df)) as.character(df[["Histology"]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}
parse_tcga_idh_from_clinical_if_available <- function(clin_raw) {
  idh_cols <- grep("idh", colnames(clin_raw), ignore.case = TRUE, value = TRUE)
  if (length(idh_cols) == 0) return(NULL)
  idh_mat <- do.call(cbind, lapply(idh_cols, function(cc) parse_idh(clin_raw[[cc]])))
  idh_vec <- apply(idh_mat, 1, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0) return(NA_real_)
    if (any(z == 1)) return(1)
    if (any(z == 0)) return(0)
    NA_real_
  })
  as.numeric(idh_vec)
}
make_tcga_idh_from_mut <- function(mut_df) {
  if (nrow(mut_df) == 0) {
    return(data.frame(
      Patient_ID = character(0),
      Mutation_Profiled = integer(0),
      IDH_from_mut = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  raw_id_col <- first_existing_col(
    mut_df,
    candidates = c("Tumor_Sample_Barcode", "sample", "Sample", "Tumor_Sample"),
    required = TRUE,
    label = "mutation sample barcode column"
  )
  gene_col <- first_existing_col(
    mut_df,
    candidates = c("gene", "Gene", "Hugo_Symbol", "SYMBOL", "symbol"),
    required = TRUE,
    label = "mutation gene symbol column"
  )
  mut_df$Patient_ID <- tcga_patient_id(mut_df[[raw_id_col]])
  all_mut_patients <- sort(unique(mut_df$Patient_ID[!is.na(mut_df$Patient_ID) & mut_df$Patient_ID != ""]))
  idh_patients <- unique(mut_df$Patient_ID[
    toupper(as.character(mut_df[[gene_col]])) %in% c("IDH1", "IDH2")
  ])
  data.frame(
    Patient_ID = all_mut_patients,
    Mutation_Profiled = 1L,
    IDH_from_mut = ifelse(all_mut_patients %in% idh_patients, 1, 0),
    stringsAsFactors = FALSE
  )
}
make_tcga_clin <- function(clin_raw, surv_raw, mut_raw, cohort_name) {
  sample_col_clin <- first_existing_col(clin_raw, c("sample", "Sample", "submitter_id", "case_submitter_id"), TRUE, "TCGA clinical sample column")
  sample_col_surv <- first_existing_col(surv_raw, c("sample", "Sample", "submitter_id", "case_submitter_id"), TRUE, "TCGA survival sample column")
  clin_raw$Sample_ID <- tcga_sample_short(clin_raw[[sample_col_clin]])
  clin_raw$Patient_ID <- tcga_patient_id(clin_raw[[sample_col_clin]])
  surv_raw$Sample_ID <- tcga_sample_short(surv_raw[[sample_col_surv]])
  surv_raw$Patient_ID <- tcga_patient_id(surv_raw[[sample_col_surv]])
  age_col <- first_existing_col(clin_raw, c("age_at_index.demographic", "age_at_diagnosis", "Age", "age"), FALSE)
  gender_col <- first_existing_col(clin_raw, c("gender.demographic", "gender", "Gender"), FALSE)
  grade_col <- first_existing_col(clin_raw, c("tumor_grade.diagnoses", "tumor_grade", "Grade", "grade"), FALSE)
  dx_col <- first_existing_col(clin_raw, c("primary_diagnosis.diagnoses", "primary_diagnosis", "diagnosis"), FALSE)
  os_time_col <- first_existing_col(surv_raw, c("OS.time", "OS_time", "time", "days_to_death", "days_to_last_follow_up"), TRUE, "OS time")
  os_status_col <- first_existing_col(surv_raw, c("OS", "OS_status", "status", "vital_status"), TRUE, "OS status")
  idh_clin <- parse_tcga_idh_from_clinical_if_available(clin_raw)
  if (is.null(idh_clin)) {
    idh_clin <- rep(NA_real_, nrow(clin_raw))
  }
  if (length(idh_clin) != nrow(clin_raw)) {
    stop("Clinical-derived IDH vector length does not match TCGA clinical rows.")
  }
  clin_raw$IDH_from_clin <- as.numeric(idh_clin)
  clin_min <- clin_raw %>%
    group_by(Sample_ID, Patient_ID) %>%
    summarise(
      Age = if (!is.null(age_col)) dplyr::first(as_num(.data[[age_col]])) else NA_real_,
      Gender = if (!is.null(gender_col)) dplyr::first(parse_gender(.data[[gender_col]])) else NA_real_,
      Grade = if (!is.null(grade_col)) dplyr::first(parse_grade(.data[[grade_col]])) else NA_real_,
      Primary_Diagnosis = if (!is.null(dx_col)) dplyr::first(as.character(.data[[dx_col]])) else NA_character_,
      IDH_from_clin = dplyr::first(.data$IDH_from_clin),
      .groups = "drop"
    )
  surv_min <- surv_raw %>%
    group_by(Sample_ID, Patient_ID) %>%
    summarise(
      OS_time = dplyr::first(as_num(.data[[os_time_col]])),
      OS_status_raw = dplyr::first(.data[[os_status_col]]),
      .groups = "drop"
    )
  surv_min$OS_status <- as_num(surv_min$OS_status_raw)
  if (all(is.na(surv_min$OS_status))) {
    tmp <- tolower(as.character(surv_min$OS_status_raw))
    surv_min$OS_status <- NA_real_
    surv_min$OS_status[grepl("dead|deceased|1", tmp)] <- 1
    surv_min$OS_status[grepl("alive|living|0", tmp)] <- 0
  }
  idh_mut <- make_tcga_idh_from_mut(mut_raw)
  out <- surv_min %>%
    left_join(clin_min, by = c("Sample_ID", "Patient_ID")) %>%
    left_join(idh_mut, by = "Patient_ID")
  # Mutation-primary TCGA IDH inference.
  # Main exposure IDH is inferred from somatic mutation calls:
  #   IDH1/IDH2 mutation present -> IDH = 1
  #   no IDH1/IDH2 mutation record -> IDH = 0
  # Clinical IDH, if available, is retained only for diagnostic comparison.
  out$Mutation_Profiled[is.na(out$Mutation_Profiled)] <- 0L
  out$IDH_from_mut[is.na(out$IDH_from_mut)] <- 0
  out$IDH <- as.numeric(out$IDH_from_mut)
  out$IDH_Discordant_Clin_vs_Mut <- ifelse(
    is.finite(out$IDH_from_clin),
    as.numeric(out$IDH_from_clin != out$IDH_from_mut),
    NA_real_
  )
  data.frame(
    Sample_ID = out$Sample_ID,
    Patient_ID = out$Patient_ID,
    Cohort = cohort_name,
    OS_time = out$OS_time,
    OS_status = out$OS_status,
    IDH = out$IDH,
    IDH_from_mut = out$IDH_from_mut,
    IDH_from_clin = out$IDH_from_clin,
    IDH_Discordant_Clin_vs_Mut = out$IDH_Discordant_Clin_vs_Mut,
    Mutation_Profiled = out$Mutation_Profiled,
    IDH_Source = "somatic_mutation_primary",
    Age = out$Age,
    Gender = out$Gender,
    Grade = out$Grade,
    Primary_Diagnosis = out$Primary_Diagnosis,
    stringsAsFactors = FALSE
  )
}
align_clin_expr <- function(clin, expr, cohort_name) {
  clin$Sample_ID <- as.character(clin$Sample_ID)
  colnames(expr) <- as.character(colnames(expr))
  common <- intersect(clin$Sample_ID, colnames(expr))
  clin2 <- clin[match(common, clin$Sample_ID), , drop = FALSE]
  expr2 <- expr[, common, drop = FALSE]
  stopifnot(all(clin2$Sample_ID == colnames(expr2)))
  message("\nAligned cohort: ", cohort_name)
  message("  Clinical rows before align: ", nrow(clin))
  message("  Expression samples before align: ", ncol(expr))
  message("  Matched samples: ", length(common))
  list(clinical = clin2, expr = expr2, cohort = cohort_name)
}
summarise_cohort <- function(obj) {
  clin <- obj$clinical
  basic_keep <- !is.na(clin$OS_time) & clin$OS_time > 0 &
    !is.na(clin$OS_status) & clin$OS_status %in% c(0, 1) &
    !is.na(clin$IDH) & clin$IDH %in% c(0, 1) &
    !is.na(clin$Age) &
    !is.na(clin$Gender)
  clin_basic <- clin[basic_keep, , drop = FALSE]
  data.frame(
    Cohort = obj$cohort,
    N_Aligned = nrow(clin),
    N_Complete_Basic = nrow(clin_basic),
    Events = sum(clin_basic$OS_status == 1, na.rm = TRUE),
    Censored = sum(clin_basic$OS_status == 0, na.rm = TRUE),
    Event_Rate = mean(clin_basic$OS_status == 1, na.rm = TRUE),
    IDH_Mutant = sum(clin_basic$IDH == 1, na.rm = TRUE),
    IDH_Wildtype = sum(clin_basic$IDH == 0, na.rm = TRUE),
    N_With_Grade = sum(!is.na(clin_basic$Grade)),
    Median_OS_Time = median(clin_basic$OS_time, na.rm = TRUE),
    N_Expression_Genes = nrow(obj$expr),
    stringsAsFactors = FALSE
  )
}
summarise_tcga_idh_diagnostics <- function(clin) {
  if (!("IDH_from_mut" %in% colnames(clin))) {
    return(data.frame())
  }
  clin_idh_available <- "IDH_from_clin" %in% colnames(clin)
  discord_available <- "IDH_Discordant_Clin_vs_Mut" %in% colnames(clin)
  mutation_profiled_available <- "Mutation_Profiled" %in% colnames(clin)
  data.frame(
    Cohort = SOURCE_COHORT,
    IDH_Main_Definition = "somatic_mutation_primary",
    N_Clinical_Rows = nrow(clin),
    N_Mutation_Profiled = if (mutation_profiled_available) sum(clin$Mutation_Profiled == 1, na.rm = TRUE) else NA_integer_,
    N_IDH_Mutant_By_Mutation = sum(clin$IDH_from_mut == 1, na.rm = TRUE),
    N_IDH_Wildtype_By_Mutation = sum(clin$IDH_from_mut == 0, na.rm = TRUE),
    N_Clinical_IDH_Available = if (clin_idh_available) sum(is.finite(clin$IDH_from_clin)) else 0L,
    N_Clinical_Mutation_Discordant = if (discord_available) sum(clin$IDH_Discordant_Clin_vs_Mut == 1, na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------------------------
# Run preparation
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Preparing final two-cohort glioma example data\n")
cat("Project directory:", PROJECT_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("============================================================\n")
cgga_clin_raw <- read_any_table(CGGA_ARRAY_CLIN)
cgga_expr <- parse_cgga_expr(CGGA_ARRAY_EXPR, TARGET_COHORT)
cgga_clin <- make_cgga_clin(cgga_clin_raw, TARGET_COHORT)
target_obj <- align_clin_expr(cgga_clin, cgga_expr, TARGET_COHORT)
tcga_clin_raw <- read_any_table(TCGA_LGG_CLIN)
tcga_surv_raw <- read_any_table(TCGA_LGG_SURV)
tcga_mut_raw <- read_any_table(TCGA_LGG_MUT)
tcga_expr_raw <- parse_tcga_expr_ensembl_to_symbol(TCGA_LGG_EXPR, SOURCE_COHORT)
tcga_expr <- choose_primary_tcga_samples_keep_sample_id(tcga_expr_raw)
tcga_clin <- make_tcga_clin(tcga_clin_raw, tcga_surv_raw, tcga_mut_raw, SOURCE_COHORT)
source_obj <- align_clin_expr(tcga_clin, tcga_expr, SOURCE_COHORT)
tcga_idh_diagnostics <- summarise_tcga_idh_diagnostics(source_obj$clinical)
common_genes <- intersect(rownames(target_obj$expr), rownames(source_obj$expr))
common_genes <- sort(common_genes)
target_obj$expr <- target_obj$expr[common_genes, , drop = FALSE]
source_obj$expr <- source_obj$expr[common_genes, , drop = FALSE]
cohort_summary <- bind_rows(
  summarise_cohort(target_obj),
  summarise_cohort(source_obj)
)
gene_overlap_summary <- data.frame(
  Target_Cohort = TARGET_COHORT,
  Source_Cohort = SOURCE_COHORT,
  Target_Genes_Before_Common = nrow(cgga_expr),
  Source_Genes_Before_Common = nrow(tcga_expr),
  Main_Common_Genes_Target_Source = length(common_genes),
  Gene_Universe_Definition = "intersection(CGGA_mRNA_array_301, TCGA_LGG)",
  stringsAsFactors = FALSE
)
target_basic_n <- cohort_summary$N_Complete_Basic[cohort_summary$Cohort == TARGET_COHORT]
source_basic_n <- cohort_summary$N_Complete_Basic[cohort_summary$Cohort == SOURCE_COHORT]
target_events <- cohort_summary$Events[cohort_summary$Cohort == TARGET_COHORT]
source_events <- cohort_summary$Events[cohort_summary$Cohort == SOURCE_COHORT]
source_target_ratio <- data.frame(
  Target_Cohort = TARGET_COHORT,
  Source_Cohort = SOURCE_COHORT,
  Target_N = target_basic_n,
  Source_N = source_basic_n,
  Source_Target_N_Ratio = source_basic_n / target_basic_n,
  Target_Events = target_events,
  Source_Events = source_events,
  Source_Target_Event_Ratio = source_events / target_events,
  stringsAsFactors = FALSE
)
expr_summary <- bind_rows(
  summarise_expr(target_obj$expr, TARGET_COHORT),
  summarise_expr(source_obj$expr, SOURCE_COHORT)
)
aligned_obj <- list(
  target = target_obj,
  source = source_obj,
  common_genes = common_genes,
  meta = list(
    target_cohort = TARGET_COHORT,
    source_cohort = SOURCE_COHORT,
    gene_universe_definition = "intersection(CGGA_mRNA_array_301, TCGA_LGG)",
    tcga_idh_definition = "somatic_mutation_primary",
    expression_transform = "no_log2",
    results_root = RESULTS_ROOT,
    created_at = as.character(Sys.time())
  )
)
saveRDS(aligned_obj, file.path(OUT_DIR, "glioma_example_aligned_object.rds"))
write.csv(cohort_summary, file.path(OUT_DIR, "glioma_example_cohort_summary.csv"), row.names = FALSE)
write.csv(gene_overlap_summary, file.path(OUT_DIR, "glioma_example_gene_overlap_summary.csv"), row.names = FALSE)
write.csv(source_target_ratio, file.path(OUT_DIR, "glioma_example_source_target_ratio.csv"), row.names = FALSE)
write.csv(expr_summary, file.path(OUT_DIR, "glioma_example_expression_summary.csv"), row.names = FALSE)
write.csv(tcga_idh_diagnostics, file.path(OUT_DIR, "glioma_example_tcga_idh_diagnostics.csv"), row.names = FALSE)
cat("\nTCGA IDH diagnostics:\n")
print(tcga_idh_diagnostics)
cat("\nPreparation complete.\n")
cat("Main common genes:", length(common_genes), "\n")
cat("Cohort summary:\n")
print(cohort_summary)
cat("\nSource-target ratio:\n")
print(source_target_ratio)
cat("\nOutput RDS:\n", file.path(OUT_DIR, "glioma_example_aligned_object.rds"), "\n", sep = "")