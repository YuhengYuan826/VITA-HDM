# ==============================================================================
# 06_summarise_simulation_results.R
# Run after file 05:
#   Rscript 06_summarise_simulation_results.R
# Summarise and plot the final cross-platform VITA-HDM simulation results.
# Publication figure specification:
#   1. Use publication-consistent notation: |T|, rho_M, and p_m.
#   2. Avoid duplicated automatic ggplot legends and guide-box warnings.
#   3. Use panel-specific custom legends and panel-specific palettes.
#   4. Give every major panel its own legend and method-mapping note.
#   5. Preserve panel-specific coordinate ranges and manuscript ticks.
#   6. Add a concise scenario header to every exported figure.
#   7. Balance the upper two-row plots and lower one-row plots visually.
#   8. Export at a full double-column width of 178 mm.
#   9. Overwrite the existing figure files.
# ==============================================================================
# ==============================================================================
# Project setup
# ==============================================================================
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(e) NULL
)
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
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}
source(file.path(SCRIPT_DIR, "00_project_setup.R"))
PROJECT_ROOT <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_ROOT)
SIMULATION_OUTPUT_DIR <- ensure_directory(PATHS$simulation_results)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(grid)
})
# ==============================================================================
# Figure export settings
# ==============================================================================
# The figure is prepared at an approximate full double-column width.
# The additional height accommodates separate panel-specific footers while
# preserving sufficient space for the lower one-row panels.
FIGURE_WIDTH_MM <- 178
FIGURE_HEIGHT_MM <- 185
FIGURE_DPI <- 600
# A generic sans-serif family preserves cross-platform portability.
FIGURE_FONT <- "sans"
# ==============================================================================
# Shared data helpers
# ==============================================================================
read_simulation_results <- function(pattern, empty_error = NULL) {
  result_files <- list.files(
    SIMULATION_OUTPUT_DIR,
    pattern = pattern,
    full.names = TRUE
  )
  if (length(result_files) == 0) {
    if (is.null(empty_error)) {
      return(NULL)
    }
    stop(empty_error, call. = FALSE)
  }
  dplyr::bind_rows(
    lapply(result_files, utils::read.csv)
  )
}
validate_required_columns <- function(data, required_columns) {
  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "The simulation results are missing the following required columns: ",
        paste(missing_columns, collapse = ", "),
        "\nPlease rerun 05_run_simulation_study.R."
      ),
      call. = FALSE
    )
  }
  invisible(data)
}
group_simulation_results <- function(data) {
  data %>%
    dplyr::group_by(
      Design,
      p_m,
      rho,
      Cens_Rate,
      Method,
      N_Source
    )
}
arrange_simulation_results <- function(data) {
  data %>%
    dplyr::arrange(
      Design,
      p_m,
      rho,
      Cens_Rate,
      Method,
      N_Source
    )
}
mean_na <- function(x) {
  valid <- !is.na(x)
  if (!any(valid)) {
    return(NA_real_)
  }
  mean(x[valid])
}
sd_na <- function(x) {
  valid <- !is.na(x)
  if (sum(valid) < 2L) {
    return(NA_real_)
  }
  stats::sd(x[valid])
}
rmse_na <- function(x) {
  valid <- !is.na(x)
  if (!any(valid)) {
    return(NA_real_)
  }
  sqrt(mean(x[valid]^2))
}
summarise_main_results <- function(data) {
  data %>%
    group_simulation_results() %>%
    dplyr::summarise(
      N_Sims = dplyr::n(),
      Actual_Censor_Mean = mean_na(Actual_Censor_Train),
      Recall_Mean = mean_na(Recall),
      FDR_Mean = mean_na(FDR),
      F1_Mean = mean_na(F1),
      Cindex_Mean = mean_na(Cindex),
      Recall_SD = sd_na(Recall),
      FDR_SD = sd_na(FDR),
      F1_SD = sd_na(F1),
      Cindex_SD = sd_na(Cindex),
      DE_SD = sd_na(DE_Est),
      IE_SD = sd_na(IE_Est),
      MP_SD = sd_na(MP_Est),
      DE_RMSE = rmse_na(DE_Error),
      IE_RMSE = rmse_na(IE_Error),
      MP_RMSE = rmse_na(MP_Error),
      Coverage_Active_Mean = mean_na(Coverage_Active),
      Coverage_Null_Mean = mean_na(Coverage_Null),
      Coverage_Active_SD = sd_na(Coverage_Active),
      Coverage_Null_SD = sd_na(Coverage_Null),
      .groups = "drop"
    ) %>%
    arrange_simulation_results()
}
summarise_sensitivity_results <- function(data) {
  data %>%
    group_simulation_results() %>%
    dplyr::summarise(
      N_Sims = dplyr::n(),
      Actual_Censor_Mean = mean_na(Actual_Censor_Train),
      Recall_Mean = mean_na(Recall),
      FDR_Mean = mean_na(FDR),
      F1_Mean = mean_na(F1),
      Cindex_Mean = mean_na(Cindex),
      DE_RMSE = rmse_na(DE_Error),
      IE_RMSE = rmse_na(IE_Error),
      MP_RMSE = rmse_na(MP_Error),
      Coverage_Active_Mean = mean_na(Coverage_Active),
      Coverage_Null_Mean = mean_na(Coverage_Null),
      .groups = "drop"
    ) %>%
    arrange_simulation_results()
}
write_summary_csv <- function(data, file_name) {
  utils::write.csv(
    data,
    file.path(SIMULATION_OUTPUT_DIR, file_name),
    row.names = FALSE
  )
}
# ==============================================================================
# Main and sensitivity summaries
# ==============================================================================
required_columns <- c(
  "Design",
  "p_m",
  "rho",
  "Cens_Rate",
  "Method",
  "N_Source",
  "Recall",
  "FDR",
  "F1",
  "Cindex",
  "DE_Error",
  "IE_Error",
  "MP_Error",
  "DE_Est",
  "IE_Est",
  "MP_Est",
  "Coverage_Active",
  "Coverage_Null",
  "Actual_Censor_Train"
)
sensitivity_required_columns <- c(
  "Design",
  "p_m",
  "rho",
  "Cens_Rate",
  "Method",
  "N_Source",
  "Recall",
  "FDR",
  "F1",
  "Cindex",
  "DE_Error",
  "IE_Error",
  "MP_Error",
  "Coverage_Active",
  "Coverage_Null",
  "Actual_Censor_Train"
)
raw <- read_simulation_results(
  pattern = "^Sim_Result_.*\\.csv$",
  empty_error = "No main simulation result files were found."
)
validate_required_columns(
  raw,
  required_columns
)
summary_df <- summarise_main_results(raw)
write_summary_csv(
  summary_df,
  "Final_Simulation_Summary_ForPlotting.csv"
)
sens_raw <- read_simulation_results(
  pattern = "^Sensitivity_StructOnly_.*\\.csv$"
)
if (!is.null(sens_raw)) {
  validate_required_columns(
    sens_raw,
    sensitivity_required_columns
  )
  sens_summary <- summarise_sensitivity_results(sens_raw)
  write_summary_csv(
    sens_summary,
    "StructureOnly_Transfer_Sensitivity_Summary.csv"
  )
}
# ==============================================================================
# Plot data, labels, and visual constants
# ==============================================================================
PM_LEVELS <- c(
  "pm1000",
  "pm2000"
)
RHO_LEVELS <- c(
  "rho0",
  "rho03",
  "rho06"
)
EFFECT_LEVELS <- c(
  "Direct effect",
  "Indirect effect",
  "Mediated proportion"
)
TOP_PLOT_ROW_WEIGHT <- 1.70
PANEL_FOOTER_ROW_WEIGHT <- 0.21
BOTTOM_PLOT_ROW_WEIGHT <- 1.00
MAIN_LAYOUT_REL_HEIGHTS <- c(
  TOP_PLOT_ROW_WEIGHT,
  PANEL_FOOTER_ROW_WEIGHT,
  BOTTOM_PLOT_ROW_WEIGHT,
  PANEL_FOOTER_ROW_WEIGHT
)
prepare_plot_data <- function(data) {
  prepared <- data %>%
    dplyr::filter(
      # Construct a single ordered curve:
      #   N_Source = 0: HDM
      #   N_Source > 0: VITA-HDM
      # This prevents multiple methods at the same x coordinate from being
      # connected when Method is not mapped to a visual aesthetic.
      (Method == "HDM" & N_Source == 0) |
        (Method == "VITA-HDM" & N_Source > 0),
      p_m %in% c(1000, 2000),
      abs(rho - 0.0) < 1e-10 |
        abs(rho - 0.3) < 1e-10 |
        abs(rho - 0.6) < 1e-10
    ) %>%
    dplyr::mutate(
      p_m_label = dplyr::case_when(
        p_m == 1000 ~ "pm1000",
        p_m == 2000 ~ "pm2000",
        TRUE ~ NA_character_
      ),
      p_m_label = factor(
        p_m_label,
        levels = PM_LEVELS
      ),
      rho_label = dplyr::case_when(
        abs(rho - 0.0) < 1e-10 ~ "rho0",
        abs(rho - 0.3) < 1e-10 ~ "rho03",
        abs(rho - 0.6) < 1e-10 ~ "rho06",
        TRUE ~ NA_character_
      ),
      rho_label = factor(
        rho_label,
        levels = RHO_LEVELS
      )
    ) %>%
    dplyr::arrange(
      Design,
      Cens_Rate,
      p_m_label,
      rho_label,
      N_Source
    )
  if (nrow(prepared) == 0) {
    stop(
      paste0(
        "No data remain after applying the plotting method rule.\n",
        "The expected rule is HDM for N_Source = 0 and ",
        "VITA-HDM for N_Source > 0."
      ),
      call. = FALSE
    )
  }
  duplicate_keys <- prepared %>%
    dplyr::count(
      Design,
      Cens_Rate,
      p_m_label,
      rho_label,
      N_Source,
      name = "n_rows"
    ) %>%
    dplyr::filter(n_rows > 1L)
  if (nrow(duplicate_keys) > 0) {
    stop(
      paste0(
        "Duplicate plotting keys were detected. Each combination of ",
        "Design, Cens_Rate, p_m, rho, and N_Source must contain ",
        "exactly one row after method filtering."
      ),
      call. = FALSE
    )
  }
  prepared
}
plot_data <- prepare_plot_data(summary_df)
pm_labeller <- ggplot2::as_labeller(
  c(
    pm1000 = "p[m] == 1000",
    pm2000 = "p[m] == 2000"
  ),
  default = ggplot2::label_parsed
)
rho_labeller <- ggplot2::as_labeller(
  c(
    rho0 = "rho[M] == 0",
    rho03 = "rho[M] == 0.3",
    rho06 = "rho[M] == 0.6"
  ),
  default = ggplot2::label_parsed
)
# Panels A, C, and D encode mediator correlation.
color_rho <- c(
  rho0 = "#3F3F3F",
  rho03 = "#0072B2",
  rho06 = "#D55E00"
)
shape_rho <- c(
  rho0 = 16,
  rho03 = 17,
  rho06 = 15
)
linetype_rho <- c(
  rho0 = "solid",
  rho03 = "dotdash",
  rho06 = "longdash"
)
# Panel B encodes effect-summary type.
color_eff <- c(
  "Direct effect" = "#111111",
  "Indirect effect" = "#CC79A7",
  "Mediated proportion" = "#009E73"
)
shape_eff <- c(
  "Direct effect" = 18,
  "Indirect effect" = 17,
  "Mediated proportion" = 16
)
linetype_eff <- c(
  "Direct effect" = "solid",
  "Indirect effect" = "longdash",
  "Mediated proportion" = "dotdash"
)
pd_rho <- ggplot2::position_dodge(width = 0.14)
pd_eff <- ggplot2::position_dodge(width = 0.20)
# Use the Unicode mathematical-script capital T through its code point.
x_lab_expr <- "|\U0001D4AF| (Number of source domains)"
method_note_text <- paste0(
  "HDM: |\U0001D4AF| = 0    ",
  "VITA-HDM: |\U0001D4AF| = 1, 2, 3"
)
panel_heading <- function(letter, title) {
  paste0(letter, "   ", title)
}
# ==============================================================================
# Publication theme and reusable plot layers
# ==============================================================================
pub_theme <- ggplot2::theme_classic(
  base_size = 9.0,
  base_family = FIGURE_FONT
) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    panel.border = ggplot2::element_rect(
      color = "#B9BDC1",
      linewidth = 0.34,
      fill = NA
    ),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(
      color = "#5C6063",
      linewidth = 0.32
    ),
    axis.ticks.length = grid::unit(1.05, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E5E7E9",
      linewidth = 0.30,
      linetype = "solid"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 9.8,
      lineheight = 0.98,
      hjust = 0,
      color = "#1D1F20",
      margin = ggplot2::margin(b = 4.2)
    ),
    axis.title = ggplot2::element_text(
      face = "plain",
      size = 8.4,
      color = "#252729"
    ),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 3.0)
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 3.0)
    ),
    axis.text = ggplot2::element_text(
      size = 7.8,
      color = "#303234"
    ),
    axis.text.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 1.5)
    ),
    axis.text.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 1.0)
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F2F3F4",
      color = "#B9BDC1",
      linewidth = 0.34
    ),
    strip.text = ggplot2::element_text(
      face = "plain",
      size = 8.0,
      lineheight = 0.95,
      color = "#252729",
      margin = ggplot2::margin(
        t = 2.4,
        r = 2.6,
        b = 2.4,
        l = 2.6
      )
    ),
    panel.spacing = grid::unit(1.20, "mm"),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      t = 3.5,
      r = 3.8,
      b = 1.8,
      l = 3.8
    )
  )
series_geoms <- function(position, line_width, point_size) {
  list(
    ggplot2::geom_line(
      linewidth = line_width,
      position = position,
      lineend = "round",
      na.rm = TRUE
    ),
    ggplot2::geom_point(
      size = point_size,
      position = position,
      stroke = 0.55,
      na.rm = TRUE
    )
  )
}
rho_scales <- function() {
  list(
    ggplot2::scale_color_manual(
      values = color_rho,
      limits = RHO_LEVELS,
      breaks = RHO_LEVELS,
      drop = FALSE
    ),
    ggplot2::scale_shape_manual(
      values = shape_rho,
      limits = RHO_LEVELS,
      breaks = RHO_LEVELS,
      drop = FALSE
    ),
    ggplot2::scale_linetype_manual(
      values = linetype_rho,
      limits = RHO_LEVELS,
      breaks = RHO_LEVELS,
      drop = FALSE
    )
  )
}
effect_scales <- function() {
  list(
    ggplot2::scale_color_manual(
      values = color_eff,
      limits = EFFECT_LEVELS,
      breaks = EFFECT_LEVELS,
      drop = FALSE
    ),
    ggplot2::scale_shape_manual(
      values = shape_eff,
      limits = EFFECT_LEVELS,
      breaks = EFFECT_LEVELS,
      drop = FALSE
    ),
    ggplot2::scale_linetype_manual(
      values = linetype_eff,
      limits = EFFECT_LEVELS,
      breaks = EFFECT_LEVELS,
      drop = FALSE
    )
  )
}
hidden_series_guides <- function() {
  ggplot2::guides(
    color = "none",
    shape = "none",
    linetype = "none"
  )
}
source_count_scale <- function() {
  ggplot2::scale_x_continuous(
    breaks = 0:3,
    minor_breaks = NULL,
    expand = ggplot2::expansion(add = 0.16)
  )
}
rho_series_layers <- function() {
  c(
    series_geoms(
      position = pd_rho,
      line_width = 0.62,
      point_size = 1.80
    ),
    rho_scales(),
    list(
      source_count_scale(),
      hidden_series_guides()
    )
  )
}
effect_series_layers <- function() {
  c(
    series_geoms(
      position = pd_eff,
      line_width = 0.64,
      point_size = 1.86
    ),
    effect_scales(),
    list(
      source_count_scale(),
      hidden_series_guides()
    )
  )
}
# ==============================================================================
# Panel functions
# ==============================================================================
draw_panel_A <- function(data) {
  long_data <- data %>%
    dplyr::select(
      N_Source,
      p_m_label,
      rho_label,
      Recall_Mean,
      FDR_Mean
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        Recall_Mean,
        FDR_Mean
      ),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Metric = factor(
        Metric,
        levels = c(
          "Recall_Mean",
          "FDR_Mean"
        ),
        labels = c(
          "Power (recall)",
          "False discovery rate"
        )
      )
    ) %>%
    dplyr::arrange(
      Metric,
      p_m_label,
      rho_label,
      N_Source
    )
  ggplot2::ggplot(
    long_data,
    ggplot2::aes(
      x = N_Source,
      y = Value,
      group = rho_label,
      color = rho_label,
      shape = rho_label,
      linetype = rho_label
    )
  ) +
    ggplot2::facet_grid(
      Metric ~ p_m_label,
      switch = "y",
      scales = "free_y",
      drop = FALSE,
      labeller = ggplot2::labeller(
        p_m_label = pm_labeller
      )
    ) +
    rho_series_layers() +
    ggplot2::labs(
      y = NULL,
      x = x_lab_expr,
      title = panel_heading(
        "A",
        "Selection performance"
      )
    ) +
    pub_theme +
    ggplot2::theme(
      strip.background.y = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y = ggplot2::element_text(
        size = 8.2,
        face = "bold",
        color = "#252729"
      )
    )
}
draw_panel_B <- function(data) {
  long_rmse <- data %>%
    dplyr::select(
      N_Source,
      p_m_label,
      rho_label,
      DE_RMSE,
      IE_RMSE,
      MP_RMSE
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        DE_RMSE,
        IE_RMSE,
        MP_RMSE
      ),
      names_to = "Type",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Type = factor(
        Type,
        levels = c(
          "DE_RMSE",
          "IE_RMSE",
          "MP_RMSE"
        ),
        labels = EFFECT_LEVELS
      )
    ) %>%
    dplyr::arrange(
      p_m_label,
      rho_label,
      Type,
      N_Source
    )
  ggplot2::ggplot(
    long_rmse,
    ggplot2::aes(
      x = N_Source,
      y = Value,
      group = Type,
      color = Type,
      shape = Type,
      linetype = Type
    )
  ) +
    ggplot2::facet_grid(
      p_m_label ~ rho_label,
      scales = "free_y",
      drop = FALSE,
      labeller = ggplot2::labeller(
        p_m_label = pm_labeller,
        rho_label = rho_labeller
      )
    ) +
    effect_series_layers() +
    ggplot2::scale_y_continuous(
      breaks = c(
        0.25,
        0.50,
        0.75,
        1.00
      ),
      minor_breaks = NULL,
      expand = ggplot2::expansion(
        mult = c(0.04, 0.08)
      )
    ) +
    ggplot2::labs(
      y = "Root mean squared error",
      x = x_lab_expr,
      title = panel_heading(
        "B",
        "Estimation accuracy"
      )
    ) +
    pub_theme
}
draw_panel_C <- function(data) {
  nominal_label <- data.frame(
    p_m_label = factor(
      "pm2000",
      levels = PM_LEVELS
    ),
    N_Source = Inf,
    Coverage_Active_Mean = 0.95,
    label = "Nominal 0.95",
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = N_Source,
      y = Coverage_Active_Mean,
      group = rho_label,
      color = rho_label,
      shape = rho_label,
      linetype = rho_label
    )
  ) +
    ggplot2::facet_grid(
      . ~ p_m_label,
      drop = FALSE,
      labeller = ggplot2::labeller(
        p_m_label = pm_labeller
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0.95,
      linetype = "dotted",
      linewidth = 0.58,
      color = "#6D7174"
    ) +
    rho_series_layers() +
    ggplot2::geom_text(
      data = nominal_label,
      mapping = ggplot2::aes(
        x = N_Source,
        y = Coverage_Active_Mean,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.06,
      vjust = -0.38,
      family = FIGURE_FONT,
      size = 2.30,
      color = "#626669"
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(
        0,
        0.25,
        0.50,
        0.75
      ),
      minor_breaks = NULL,
      expand = ggplot2::expansion(
        mult = c(0, 0.01)
      )
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, 1.05),
      clip = "on"
    ) +
    ggplot2::labs(
      y = "Active coverage",
      x = x_lab_expr,
      title = panel_heading(
        "C",
        "Variational interval coverage"
      )
    ) +
    pub_theme
}
draw_panel_D <- function(data) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = N_Source,
      y = Cindex_Mean,
      group = rho_label,
      color = rho_label,
      shape = rho_label,
      linetype = rho_label
    )
  ) +
    ggplot2::facet_grid(
      . ~ p_m_label,
      drop = FALSE,
      labeller = ggplot2::labeller(
        p_m_label = pm_labeller
      )
    ) +
    rho_series_layers() +
    ggplot2::scale_y_continuous(
      breaks = seq(
        0.86,
        0.90,
        by = 0.01
      ),
      labels = function(x) {
        sprintf("%.2f", x)
      },
      minor_breaks = NULL,
      expand = ggplot2::expansion(
        mult = c(0.05, 0.08)
      )
    ) +
    ggplot2::labs(
      y = "Concordance index",
      x = x_lab_expr,
      title = panel_heading(
        "D",
        "Predictive performance"
      )
    ) +
    pub_theme
}
# ==============================================================================
# Panel-specific footer plots
# ==============================================================================
make_panel_footer <- function(keys,
                              x,
                              labels,
                              scale_layers,
                              segment_start,
                              segment_end,
                              point_offset,
                              text_offset,
                              point_size,
                              line_width,
                              parse_labels = FALSE) {
  legend_data <- data.frame(
    x = x,
    y = rep(0.70, length(x)),
    legend_key = factor(
      keys,
      levels = keys
    ),
    label = labels,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(legend_data) +
    ggplot2::annotate(
      "segment",
      x = 0.01,
      xend = 0.99,
      y = 0.98,
      yend = 0.98,
      color = "#D7D9DB",
      linewidth = 0.28
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = x - segment_start,
        xend = x + segment_end,
        y = y,
        yend = y,
        color = legend_key,
        linetype = legend_key
      ),
      linewidth = line_width,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = x - point_offset,
        y = y,
        color = legend_key,
        shape = legend_key
      ),
      size = point_size,
      stroke = 0.55
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = x + text_offset,
        y = y,
        label = label
      ),
      parse = parse_labels,
      hjust = 0,
      vjust = 0.5,
      family = FIGURE_FONT,
      color = "#303234",
      size = 2.42
    ) +
    ggplot2::annotate(
      "text",
      x = 0.50,
      y = 0.20,
      label = method_note_text,
      hjust = 0.5,
      vjust = 0.5,
      family = FIGURE_FONT,
      fontface = "plain",
      color = "#4D5052",
      size = 2.20
    ) +
    scale_layers +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      clip = "on"
    ) +
    hidden_series_guides() +
    ggplot2::theme_void(
      base_size = 9,
      base_family = FIGURE_FONT
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      plot.margin = ggplot2::margin(
        t = 0,
        r = 3,
        b = 0,
        l = 3
      )
    )
}
make_rho_panel_footer <- function() {
  make_panel_footer(
    keys = RHO_LEVELS,
    x = c(
      0.14,
      0.45,
      0.76
    ),
    labels = c(
      "rho[M] == 0",
      "rho[M] == 0.3",
      "rho[M] == 0.6"
    ),
    scale_layers = rho_scales(),
    segment_start = 0.045,
    segment_end = 0.014,
    point_offset = 0.016,
    text_offset = 0.026,
    point_size = 1.72,
    line_width = 0.60,
    parse_labels = TRUE
  )
}
make_effect_panel_footer <- function() {
  make_panel_footer(
    keys = EFFECT_LEVELS,
    x = c(
      0.13,
      0.37,
      0.67
    ),
    labels = c(
      "Direct",
      "Indirect",
      "Mediated proportion"
    ),
    scale_layers = effect_scales(),
    segment_start = 0.045,
    segment_end = 0.014,
    point_offset = 0.016,
    text_offset = 0.028,
    point_size = 1.78,
    line_width = 0.62,
    parse_labels = FALSE
  )
}
make_top_footer_row <- function() {
  cowplot::plot_grid(
    # Panel A footer
    make_rho_panel_footer(),
    # Panel B footer
    make_effect_panel_footer(),
    ncol = 2,
    rel_widths = c(1, 1),
    align = "h",
    axis = "tb"
  )
}
make_bottom_footer_row <- function() {
  cowplot::plot_grid(
    # Panel C footer
    make_rho_panel_footer(),
    # Panel D footer
    make_rho_panel_footer(),
    ncol = 2,
    rel_widths = c(1, 1),
    align = "h",
    axis = "tb"
  )
}
# ==============================================================================
# Scenario header
# ==============================================================================
make_scenario_header <- function(design, cens_rate) {
  header_text <- paste0(
    design,
    " design  \u00B7  ",
    sprintf(
      "%.0f%% censoring",
      100 * cens_rate
    )
  )
  cowplot::ggdraw() +
    cowplot::draw_grob(
      grid::rectGrob(
        gp = grid::gpar(
          fill = "#F3F4F5",
          col = "#D5D8DA",
          lwd = 0.65
        )
      )
    ) +
    cowplot::draw_label(
      header_text,
      x = 0.014,
      y = 0.50,
      hjust = 0,
      vjust = 0.5,
      fontfamily = FIGURE_FONT,
      fontface = "bold",
      size = 9.0,
      color = "#252729"
    )
}
# ==============================================================================
# Composite figure builder and export
# ==============================================================================
make_simulation_figure <- function(data,
                                   design,
                                   cens_rate,
                                   output_base,
                                   output_dir = ".") {
  scenario_data <- data %>%
    dplyr::filter(
      Design == design,
      abs(Cens_Rate - cens_rate) < 1e-10
    ) %>%
    dplyr::arrange(
      p_m_label,
      rho_label,
      N_Source
    )
  if (nrow(scenario_data) == 0) {
    stop(
      paste0(
        "No data were found for Design = ",
        design,
        " and Cens_Rate = ",
        cens_rate,
        "."
      ),
      call. = FALSE
    )
  }
  panel_A <- draw_panel_A(scenario_data)
  panel_B <- draw_panel_B(scenario_data)
  panel_C <- draw_panel_C(scenario_data)
  panel_D <- draw_panel_D(scenario_data)
  aligned_panels <- cowplot::align_plots(
    panel_A,
    panel_B,
    panel_C,
    panel_D,
    align = "hv",
    axis = "tblr"
  )
  top_plot_row <- cowplot::plot_grid(
    aligned_panels[[1]],
    aligned_panels[[2]],
    ncol = 2,
    rel_widths = c(1, 1)
  )
  top_footer_row <- make_top_footer_row()
  bottom_plot_row <- cowplot::plot_grid(
    aligned_panels[[3]],
    aligned_panels[[4]],
    ncol = 2,
    rel_widths = c(1, 1)
  )
  bottom_footer_row <- make_bottom_footer_row()
  main_grid <- cowplot::plot_grid(
    top_plot_row,
    top_footer_row,
    bottom_plot_row,
    bottom_footer_row,
    ncol = 1,
    rel_heights = MAIN_LAYOUT_REL_HEIGHTS
  )
  scenario_header <- make_scenario_header(
    design = design,
    cens_rate = cens_rate
  )
  final_plot <- cowplot::plot_grid(
    scenario_header,
    main_grid,
    ncol = 1,
    rel_heights = c(
      0.050,
      1.000
    )
  )
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  # These paths intentionally use the original filenames.
  # ggsave() overwrites existing files automatically.
  png_path <- file.path(
    output_dir,
    paste0(output_base, ".png")
  )
  pdf_path <- file.path(
    output_dir,
    paste0(output_base, ".pdf")
  )
  # Prefer ragg when available, while retaining a standard PNG fallback.
  png_device <- if (
    requireNamespace("ragg", quietly = TRUE)
  ) {
    ragg::agg_png
  } else {
    grDevices::png
  }
  ggplot2::ggsave(
    filename = png_path,
    plot = final_plot,
    width = FIGURE_WIDTH_MM,
    height = FIGURE_HEIGHT_MM,
    units = "mm",
    dpi = FIGURE_DPI,
    device = png_device,
    bg = "white",
    limitsize = FALSE
  )
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  ggplot2::ggsave(
    filename = pdf_path,
    plot = final_plot,
    width = FIGURE_WIDTH_MM,
    height = FIGURE_HEIGHT_MM,
    units = "mm",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )
  message("Saved: ", png_path)
  message("Saved: ", pdf_path)
  invisible(final_plot)
}
# ==============================================================================
# Export all four final simulation figures
# ==============================================================================
figure_scenarios <- data.frame(
  design = c(
    "Homogeneous",
    "Homogeneous",
    "Heterogeneous",
    "Heterogeneous"
  ),
  cens_rate = c(
    0.20,
    0.40,
    0.20,
    0.40
  ),
  output_base = c(
    "Figure_FinalBalanced_Homogeneous_Cens20",
    "Figure_FinalBalanced_Homogeneous_Cens40",
    "Figure_FinalBalanced_Heterogeneous_Cens20",
    "Figure_FinalBalanced_Heterogeneous_Cens40"
  ),
  stringsAsFactors = FALSE
)
for (scenario_index in seq_len(nrow(figure_scenarios))) {
  scenario <- figure_scenarios[scenario_index, ]
  make_simulation_figure(
    data = plot_data,
    design = scenario$design,
    cens_rate = scenario$cens_rate,
    output_base = scenario$output_base,
    output_dir = SIMULATION_OUTPUT_DIR
  )
}
message("Simulation summarisation and figure export completed successfully.")