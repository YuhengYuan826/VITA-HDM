# ==============================================================================
# 00_project_setup.R
# Shared repository paths and dependency inventory.
# Scientific parameters remain in the numbered analysis files.
# ==============================================================================
find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "00_project_setup.R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop(
    "Cannot locate 00_project_setup.R from: ",
    start,
    call. = FALSE
  )
}
project_paths <- function(root = find_project_root()) {
  data_root <- Sys.getenv(
    "VITA_HDM_DATA_DIR",
    unset = file.path(root, "data")
  )
  results_root <- Sys.getenv(
    "VITA_HDM_RESULTS_DIR",
    unset = file.path(root, "results")
  )
  data_root <- normalizePath(
    data_root,
    winslash = "/",
    mustWork = FALSE
  )
  results_root <- normalizePath(
    results_root,
    winslash = "/",
    mustWork = FALSE
  )
  case_root <- file.path(results_root, "case_study")
  sensitivity_root <- file.path(case_root, "sensitivity")
  list(
    root = root,
    data = data_root,
    cgga_data = file.path(data_root, "CGGA"),
    tcga_data = file.path(data_root, "TCGA"),
    results = results_root,
    simulation_results = file.path(results_root, "simulation"),
    case_results = case_root,
    case_data = file.path(case_root, "data"),
    discovery = file.path(case_root, "discovery"),
    fixed_signature_cv = file.path(case_root, "fixed_signature_cv"),
    conditional_effects = file.path(case_root, "conditional_effects"),
    case_summary = file.path(case_root, "summary"),
    sensitivity = sensitivity_root,
    sensitivity_tables = file.path(sensitivity_root, "tables"),
    sensitivity_figures = file.path(sensitivity_root, "figures"),
    sensitivity_private = file.path(
      sensitivity_root,
      "private_participant_level"
    ),
    sensitivity_final = file.path(sensitivity_root, "final_outputs")
  )
}
ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Cannot create output directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
required_cran_packages <- c(
  "MASS",
  "survival",
  "dplyr",
  "tidyr",
  "ggplot2",
  "cowplot",
  "doParallel",
  "foreach",
  "data.table",
  "stringr",
  "scales",
  "patchwork"
)
required_bioconductor_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)
optional_packages <- c(
  "ggrepel"
)
dependency_status <- function() {
  package_names <- c(
    required_cran_packages,
    required_bioconductor_packages,
    optional_packages
  )
  package_type <- c(
    rep("CRAN", length(required_cran_packages)),
    rep("Bioconductor", length(required_bioconductor_packages)),
    rep("Optional", length(optional_packages))
  )
  required <- package_type != "Optional"
  installed <- vapply(
    package_names,
    function(package_name) {
      requireNamespace(package_name, quietly = TRUE)
    },
    logical(1)
  )
  version <- vapply(
    package_names,
    function(package_name) {
      if (!requireNamespace(package_name, quietly = TRUE)) {
        return(NA_character_)
      }
      as.character(utils::packageVersion(package_name))
    },
    character(1)
  )
  data.frame(
    Package = package_names,
    Type = package_type,
    Required = required,
    Installed = installed,
    Version = version,
    stringsAsFactors = FALSE
  )
}
assert_required_dependencies <- function() {
  status <- dependency_status()
  missing <- status$Package[status$Required & !status$Installed]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(status)
}