# ==============================================================================
# 90_check_environment.R
# Optional preflight: Rscript 90_check_environment.R
# Check and record the R environment without running an analysis.
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

metadata_dir <- ensure_directory(file.path(PATHS$results, "run_metadata"))
status <- dependency_status()

utils::write.csv(
  status,
  file.path(metadata_dir, "package_versions.csv"),
  row.names = FALSE,
  na = ""
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(metadata_dir, "session_info.txt")
)
utils::write.csv(
  data.frame(
    Capability = c("cairo", "png", "jpeg", "tiff"),
    Available = vapply(
      c("cairo", "png", "jpeg", "tiff"),
      capabilities,
      FUN.VALUE = logical(1)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(metadata_dir, "graphics_capabilities.csv"),
  row.names = FALSE
)

print(status)
assert_required_dependencies()

message("Environment metadata written to: ", metadata_dir)
