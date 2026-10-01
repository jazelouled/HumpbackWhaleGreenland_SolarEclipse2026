###############################################################################
# HUMPBACK WHALE SOLAR ECLIPSE — FINAL REPRODUCIBLE ANALYSIS
#
# One script to run from the raw Wildlife Computers exports. It creates a new
# `outputs_final` directory containing all data products, models, Tables S1–S8,
# supplementary figures, a report, and manuscript Figure 4.
#
# Statistical framework
#   * Candidate residual-correlation structures are compared by ML separately
#     for dive duration and maximum dive depth.
#   * The selected structure for each response is then refitted by REML.
#   * The selected duration structure is reused for duration-model sensitivities
#     and episode analyses; the selected depth structure is reused for depth
#     models and the baseline depth plot.
###############################################################################

# =============================================================================
# 0. USER SETTINGS
# =============================================================================

# Change only this line if the project is stored elsewhere.
project_dir <- "/Users/jazelouled-cheikhbonan/Dropbox/2026_Humpback_Eclipse_RdeStephanis"

input_dir <- file.path(project_dir, "00input", "260680")
output_dir <- file.path(project_dir, "outputs_final")

behaviour_file <- file.path(input_dir, "260680-Behavior.csv")
locations_file <- file.path(input_dir, "260680-Locations.csv")

# Set to FALSE only for a quick development run. The submitted analysis should
# retain the moving-block bootstrap with 10,000 draws.
run_block_bootstrap <- TRUE
bootstrap_replicates <- 10000

# =============================================================================
# 1. PACKAGES AND OUTPUT FOLDERS
# =============================================================================

required_packages <- c(
  "mgcv", "nlme", "suncalc", "dplyr", "ggplot2", "grid",
  "gridExtra", "svglite", "scales"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))

if (!file.exists(behaviour_file) || !file.exists(locations_file)) {
  stop(
    "Raw input files were not found. Check `project_dir` and verify that\n",
    behaviour_file, "\nand\n", locations_file, "\nexist."
  )
}

data_dir <- file.path(output_dir, "data")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
audit_dir <- file.path(output_dir, "audit")

for (directory in c(output_dir, data_dir, table_dir, figure_dir, model_dir, audit_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

write_table <- function(x, filename, directory = table_dir) {
  write.csv(x, file.path(directory, filename), row.names = FALSE, na = "")
}

# =============================================================================
# 2. STUDY CONSTANTS
# =============================================================================

deployment_date <- as.Date("2026-08-09")
last_incomplete_date <- as.Date("2026-08-20")

eclipse_start <- as.POSIXct("2026-08-12 16:37:10", tz = "UTC")
eclipse_max <- as.POSIXct("2026-08-12 17:41:00", tz = "UTC")
eclipse_end <- as.POSIXct("2026-08-12 18:42:30", tz = "UTC")
eclipse_duration_min <- as.numeric(difftime(eclipse_end, eclipse_start, units = "mins"))

principal_k <- 8
solar_margins <- c(0, 2.5, 5, 7.5, 10, Inf)
k_values <- c(6, 8, 10, 12)
window_minimum_dives <- c(3, 5, 7, 9)
bootstrap_block_dives <- c(2, 3, 4, 5, 6, 8, 10, 13)

# Local eclipse geometry for the exploratory obscuration analysis.
maximum_obscuration <- 0.94
obscuration_time_step_seconds <- 30

# =============================================================================
# 3. HELPERS
# =============================================================================

parse_wc_datetime <- function(x) {
  as.POSIXct(x, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")
}

make_correlation <- function(type) {
  switch(
    type,
    NONE = NULL,
    CAR1 = nlme::corCAR1(form = ~ global_time_hours),
    AR1 = nlme::corAR1(form = ~ global_dive_index),
    AR2 = nlme::corARMA(p = 2, q = 0, form = ~ global_dive_index),
    ARMA11 = nlme::corARMA(p = 1, q = 1, form = ~ global_dive_index),
    stop("Unknown correlation type: ", type)
  )
}

temporal_labels <- c(
  NONE = "No temporal correlation",
  CAR1 = "Continuous-time AR(1)",
  AR1 = "Discrete AR(1)",
  AR2 = "Discrete AR(2)",
  ARMA11 = "Discrete ARMA(1,1)"
)

fit_gamm <- function(formula, data, correlation_type, method = c("ML", "REML")) {
  method <- match.arg(method)
  correlation <- make_correlation(correlation_type)
  args <- list(formula = formula, data = data, method = method)
  if (!is.null(correlation)) args$correlation <- correlation
  do.call(mgcv::gamm, args)
}

extract_eclipse_effect <- function(model, response, unit) {
  coefficient_table <- summary(model$gam)$p.table
  if (!"eclipse_feclipse" %in% rownames(coefficient_table)) {
    stop("Could not find the eclipse coefficient in the fitted model.")
  }
  estimate <- coefficient_table["eclipse_feclipse", "Estimate"]
  se <- coefficient_table["eclipse_feclipse", "Std. Error"]
  data.frame(
    response = response,
    unit = unit,
    estimate = unname(estimate),
    SE = unname(se),
    lower95 = unname(estimate - 1.96 * se),
    upper95 = unname(estimate + 1.96 * se),
    p = unname(coefficient_table["eclipse_feclipse", ncol(coefficient_table)]),
    stringsAsFactors = FALSE
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p)))
}

model_aic_row <- function(model, response, correlation_type) {
  ll <- logLik(model$lme)
  data.frame(
    response = response,
    correlation_type = correlation_type,
    temporal_structure = unname(temporal_labels[correlation_type]),
    log_likelihood = as.numeric(ll),
    n_parameters = attr(ll, "df"),
    AIC = AIC(model$lme),
    stringsAsFactors = FALSE
  )
}

compare_temporal_structures <- function(formula, data, response) {
  candidate_types <- c("NONE", "CAR1", "AR1", "AR2", "ARMA11")
  models <- setNames(
    lapply(candidate_types, function(x) fit_gamm(formula, data, x, method = "ML")),
    candidate_types
  )
  comparison <- bind_rows(lapply(names(models), function(x) {
    model_aic_row(models[[x]], response, x)
  })) %>%
    arrange(AIC) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      selected = delta_AIC == 0
    )
  list(models = models, comparison = comparison, selected_type = comparison$correlation_type[1])
}

save_figure <- function(plot_or_grob, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot_or_grob,
         width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot_or_grob,
         width = width, height = height, device = cairo_pdf, bg = "white")
}

# =============================================================================
# 4. RAW DATA PROCESSING
# =============================================================================

behaviour_raw <- read.csv(behaviour_file, stringsAsFactors = FALSE)
locations_raw <- read.csv(locations_file, stringsAsFactors = FALSE)

raw_dive_records <- behaviour_raw %>%
  filter(What == "Dive") %>%
  mutate(
    start_time = parse_wc_datetime(Start),
    end_time = parse_wc_datetime(End),
    duration_from_times_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
    duration_min = duration_from_times_sec / 60,
    duration_bin_mid_sec = (DurationMin + DurationMax) / 2,
    duration_bin_width_sec = DurationMax - DurationMin,
    dive_depth_m = (DepthMin + DepthMax) / 2,
    depth_bin_width_m = DepthMax - DepthMin,
    depth_bin_relative_width_pct = 100 * depth_bin_width_m / dive_depth_m
  ) %>%
  filter(
    !is.na(start_time), !is.na(end_time),
    is.finite(duration_min), is.finite(DepthMin), is.finite(DepthMax),
    as.Date(start_time, tz = "UTC") >= deployment_date
  ) %>%
  arrange(start_time)

duration_consistency <- raw_dive_records %>%
  mutate(within_transmitted_interval = duration_from_times_sec >= DurationMin &
           duration_from_times_sec <= DurationMax) %>%
  summarise(
    n_dives = n(),
    n_within_transmitted_interval = sum(within_transmitted_interval, na.rm = TRUE),
    proportion_within_transmitted_interval = mean(within_transmitted_interval, na.rm = TRUE),
    median_absolute_difference_from_duration_midpoint_sec = median(
      abs(duration_from_times_sec - duration_bin_mid_sec), na.rm = TRUE
    )
  )

depth_resolution <- raw_dive_records %>%
  summarise(
    median_depth_bin_width_m = median(depth_bin_width_m, na.rm = TRUE),
    q95_depth_bin_width_m = quantile(depth_bin_width_m, 0.95, na.rm = TRUE),
    median_relative_depth_bin_width_pct = median(depth_bin_relative_width_pct, na.rm = TRUE),
    q95_relative_depth_bin_width_pct = quantile(depth_bin_relative_width_pct, 0.95, na.rm = TRUE)
  )

locations_good <- locations_raw %>%
  mutate(
    datetime = parse_wc_datetime(Date),
    date = as.Date(datetime, tz = "UTC"),
    Quality = as.character(Quality)
  ) %>%
  filter(
    !is.na(datetime), is.finite(Latitude), is.finite(Longitude),
    Quality %in% c("1", "2", "3")
  )

daily_position <- locations_good %>%
  group_by(date) %>%
  summarise(
    daily_lat = median(Latitude, na.rm = TRUE),
    daily_lon = median(Longitude, na.rm = TRUE),
    n_argos_good = n(),
    .groups = "drop"
  )

dives <- raw_dive_records %>%
  mutate(
    midpoint = start_time + (end_time - start_time) / 2,
    date = as.Date(midpoint, tz = "UTC")
  ) %>%
  left_join(daily_position, by = "date") %>%
  filter(is.finite(daily_lat), is.finite(daily_lon)) %>%
  arrange(midpoint)

solar_input <- data.frame(date = dives$midpoint, lat = dives$daily_lat, lon = dives$daily_lon)
solar_position <- suncalc::getSunlightPosition(data = solar_input)
dives$solar_elev_deg <- solar_position$altitude * 180 / pi
dives$solar_azimuth_deg <- solar_position$azimuth * 180 / pi

dives <- dives %>%
  mutate(
    global_dive_index = row_number(),
    global_time_hours = as.numeric(difftime(midpoint, min(midpoint), units = "hours")),
    real_gap_min = as.numeric(difftime(midpoint, lag(midpoint), units = "mins")),
    prev_duration_min = lag(duration_min),
    prev_depth_m = lag(dive_depth_m),
    eclipse_midpoint = midpoint >= eclipse_start & midpoint <= eclipse_end,
    eclipse_overlap = start_time <= eclipse_end & end_time >= eclipse_start,
    eclipse = eclipse_midpoint,
    eclipse_f = factor(ifelse(eclipse, "eclipse", "control"),
                       levels = c("control", "eclipse"))
  )

eclipse_definition_summary <- data.frame(
  definition = c("midpoint", "any_overlap"),
  n_eclipse_dives = c(sum(dives$eclipse_midpoint), sum(dives$eclipse_overlap))
)

eclipse_solar_min <- min(dives$solar_elev_deg[dives$eclipse_midpoint], na.rm = TRUE)
eclipse_solar_max <- max(dives$solar_elev_deg[dives$eclipse_midpoint], na.rm = TRUE)

build_analysis_dataset <- function(solar_margin = 0, exclude_partial_days = TRUE,
                                   depth_representation = c("midpoint", "min", "max"),
                                   eclipse_definition = c("midpoint", "overlap")) {
  depth_representation <- match.arg(depth_representation)
  eclipse_definition <- match.arg(eclipse_definition)
  dat <- dives

  dat$dive_depth_model <- switch(
    depth_representation,
    midpoint = (dat$DepthMin + dat$DepthMax) / 2,
    min = dat$DepthMin,
    max = dat$DepthMax
  )
  dat$eclipse_model <- if (eclipse_definition == "midpoint") {
    dat$eclipse_midpoint
  } else {
    dat$eclipse_overlap
  }
  dat$eclipse_f <- factor(ifelse(dat$eclipse_model, "eclipse", "control"),
                          levels = c("control", "eclipse"))

  event_solar_min <- min(dat$solar_elev_deg[dat$eclipse_model], na.rm = TRUE)
  event_solar_max <- max(dat$solar_elev_deg[dat$eclipse_model], na.rm = TRUE)

  if (exclude_partial_days) {
    dat <- dat %>% filter(date != deployment_date, date != last_incomplete_date)
  }
  if (is.finite(solar_margin)) {
    dat <- dat %>% filter(
      solar_elev_deg >= event_solar_min - solar_margin,
      solar_elev_deg <= event_solar_max + solar_margin
    )
  }
  dat %>%
    arrange(midpoint) %>%
    mutate(
      analysis_row = row_number(),
      complete_record_lag = global_dive_index - lag(global_dive_index),
      analysis_gap_min = as.numeric(difftime(midpoint, lag(midpoint), units = "mins"))
    )
}

analysis <- build_analysis_dataset()
window_analysis <- build_analysis_dataset(solar_margin = Inf)

if (sum(analysis$eclipse_model) == 0 || sum(!analysis$eclipse_model) == 0) {
  stop("The principal analysis dataset contains no eclipse or no control dives.")
}

# =============================================================================
# 5. PRINCIPAL MODELS AND TEMPORAL-STRUCTURE SELECTION
# =============================================================================

duration_formula <- duration_min ~ eclipse_f +
  s(solar_elev_deg, k = principal_k, bs = "tp") +
  s(dive_depth_model, k = principal_k, bs = "tp")

depth_formula <- dive_depth_model ~ eclipse_f +
  s(solar_elev_deg, k = principal_k, bs = "tp")

# AIC comparisons use maximum likelihood. These are the only models used to
# choose correlation type; all reported coefficients come from REML refits.
duration_temporal <- compare_temporal_structures(duration_formula, analysis, "Dive duration")
depth_temporal <- compare_temporal_structures(depth_formula, analysis, "Dive depth")

duration_correlation_type <- duration_temporal$selected_type
depth_correlation_type <- depth_temporal$selected_type

m_duration <- fit_gamm(duration_formula, analysis, duration_correlation_type, method = "REML")
m_depth <- fit_gamm(depth_formula, analysis, depth_correlation_type, method = "REML")

# Baseline solar models use the response-specific structure selected above.
m_base_duration <- fit_gamm(
  duration_min ~ s(solar_elev_deg, k = 10, bs = "tp"),
  dives, duration_correlation_type, method = "REML"
)
m_base_depth <- fit_gamm(
  dive_depth_m ~ s(solar_elev_deg, k = 10, bs = "tp"),
  dives, depth_correlation_type, method = "REML"
)

principal_duration_effect <- extract_eclipse_effect(m_duration, "Dive duration", "min")
principal_depth_effect <- extract_eclipse_effect(m_depth, "Dive depth", "m")

temporal_table <- bind_rows(duration_temporal$comparison, depth_temporal$comparison)

# =============================================================================
# 6. COMMON SUPPORT AND DURATION-MODEL SENSITIVITIES
# =============================================================================

support_eclipse <- analysis %>% filter(eclipse_model)
support_control <- analysis %>% filter(!eclipse_model)
nearest_control_distance <- vapply(support_eclipse$solar_elev_deg, function(x) {
  min(abs(support_control$solar_elev_deg - x), na.rm = TRUE)
}, numeric(1))

common_support_summary <- data.frame(
  principal_N = nrow(analysis),
  eclipse_N = nrow(support_eclipse),
  control_N = nrow(support_control),
  eclipse_solar_min = min(support_eclipse$solar_elev_deg),
  eclipse_solar_max = max(support_eclipse$solar_elev_deg),
  control_solar_min = min(support_control$solar_elev_deg),
  control_solar_max = max(support_control$solar_elev_deg),
  nearest_control_median_deg = median(nearest_control_distance),
  nearest_control_max_deg = max(nearest_control_distance)
)

nearest_control_table <- data.frame(
  eclipse_solar_elev_deg = support_eclipse$solar_elev_deg,
  nearest_control_distance_deg = nearest_control_distance
)

local_support_table <- bind_rows(lapply(c(0.25, 0.5, 1, 2), function(tolerance) {
  counts <- vapply(support_eclipse$solar_elev_deg, function(x) {
    sum(abs(support_control$solar_elev_deg - x) <= tolerance)
  }, numeric(1))
  data.frame(
    tolerance_deg = tolerance,
    minimum_controls = min(counts),
    median_controls = median(counts),
    maximum_controls = max(counts)
  )
}))

fit_duration_model <- function(data, k = principal_k, include_depth = TRUE,
                               correlation_type = duration_correlation_type) {
  formula <- if (include_depth) {
    duration_formula <- duration_min ~ eclipse_f +
      s(solar_elev_deg, k = k, bs = "tp") +
      s(dive_depth_model, k = k, bs = "tp")
  } else {
    duration_formula <- duration_min ~ eclipse_f +
      s(solar_elev_deg, k = k, bs = "tp")
  }
  fit_gamm(duration_formula, data, correlation_type, method = "REML")
}

summarise_duration_model <- function(model, analysis_name, n_dives) {
  extract_eclipse_effect(model, "Dive duration", "min") %>%
    transmute(
      analysis = analysis_name,
      n_dives = n_dives,
      estimate_min = estimate,
      lower95, upper95, p
    )
}

sensitivity_rows <- list(
  summarise_duration_model(m_duration, "Principal model", nrow(analysis))
)

# Eclipse definition and depth representation.
for (definition in c("overlap")) {
  dat <- build_analysis_dataset(eclipse_definition = definition)
  model <- fit_duration_model(dat)
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
    model, "Any eclipse overlap", nrow(dat)
  )
}

for (representation in c("min", "max")) {
  dat <- build_analysis_dataset(depth_representation = representation)
  model <- fit_duration_model(dat)
  label <- if (representation == "min") "Depth minimum" else "Depth maximum"
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(model, label, nrow(dat))
}

# Solar support, including the exact range as a transparent duplicate of the
# principal analysis and wider control definitions.
for (margin in solar_margins) {
  dat <- build_analysis_dataset(solar_margin = margin)
  model <- fit_duration_model(dat)
  label <- if (margin == 0) {
    "Exact eclipse solar range"
  } else if (is.infinite(margin)) {
    "No solar restriction"
  } else {
    paste0("Solar range plus or minus ", margin, " degrees")
  }
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(model, label, nrow(dat))
}

# Partial days, depth adjustment, residual-correlation choice and smooth basis.
partial_day_data <- build_analysis_dataset(exclude_partial_days = FALSE)
sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
  fit_duration_model(partial_day_data), "Include partial recording days", nrow(partial_day_data)
)

without_depth <- fit_duration_model(analysis, include_depth = FALSE)
sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
  without_depth, "Solar elevation only (no depth adjustment)", nrow(analysis)
)

for (correlation_type in c("NONE", "CAR1", "AR1", "AR2", "ARMA11")) {
  model <- fit_duration_model(analysis, correlation_type = correlation_type)
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
    model, unname(temporal_labels[correlation_type]), nrow(analysis)
  )
}

for (k_value in k_values) {
  model <- fit_duration_model(analysis, k = k_value)
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
    model, paste0("Smooth basis k = ", k_value), nrow(analysis)
  )
}

# Previous-dive diagnostics retain the selected duration structure.
m_prev_duration <- fit_gamm(
  duration_min ~ eclipse_f + s(solar_elev_deg, k = principal_k, bs = "tp") +
    s(dive_depth_model, k = principal_k, bs = "tp") + prev_duration_min,
  analysis, duration_correlation_type, method = "REML"
)
m_prev_duration_depth <- fit_gamm(
  duration_min ~ eclipse_f + s(solar_elev_deg, k = principal_k, bs = "tp") +
    s(dive_depth_model, k = principal_k, bs = "tp") + prev_duration_min + prev_depth_m,
  analysis, duration_correlation_type, method = "REML"
)
sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
  m_prev_duration, "Previous dive duration", nrow(analysis)
)
sensitivity_rows[[length(sensitivity_rows) + 1]] <- summarise_duration_model(
  m_prev_duration_depth, "Previous dive duration and depth", nrow(analysis)
)

duration_sensitivity_table <- bind_rows(sensitivity_rows)

# =============================================================================
# 7. COUNTERFACTUAL DURATION EFFECT
# =============================================================================

actual_eclipse <- analysis %>% filter(eclipse_model)
counterfactual_control <- actual_eclipse
counterfactual_control$eclipse_f <- factor("control", levels = c("control", "eclipse"))
counterfactual_eclipse <- actual_eclipse
counterfactual_eclipse$eclipse_f <- factor("eclipse", levels = c("control", "eclipse"))

predicted_control <- predict(m_duration$gam, newdata = counterfactual_control, type = "response")
predicted_eclipse <- predict(m_duration$gam, newdata = counterfactual_eclipse, type = "response")

counterfactual_table <- data.frame(
  mean_predicted_non_eclipse_min = mean(predicted_control),
  mean_predicted_eclipse_min = mean(predicted_eclipse),
  mean_observed_eclipse_min = mean(actual_eclipse$duration_min),
  adjusted_difference_min = mean(predicted_eclipse) - mean(predicted_control),
  reduction_percent = 100 * abs(mean(predicted_eclipse) - mean(predicted_control)) / mean(predicted_control)
)

# =============================================================================
# 8. 125-MIN DURATION-EPISODE ANALYSIS
# =============================================================================

m_window_reference <- fit_gamm(
  duration_min ~ s(solar_elev_deg, k = principal_k, bs = "tp") +
    s(dive_depth_model, k = principal_k, bs = "tp"),
  window_analysis, duration_correlation_type, method = "REML"
)

window_analysis$expected_duration <- predict(
  m_window_reference$gam, newdata = window_analysis, type = "response"
)
window_analysis$adjusted_duration_residual <-
  window_analysis$duration_min - window_analysis$expected_duration

actual_eclipse_window <- window_analysis %>% filter(eclipse_model)
actual_mean_residual <- mean(actual_eclipse_window$adjusted_duration_residual)
eclipse_n_dives <- nrow(actual_eclipse_window)

build_windows <- function(data, minimum_dives = 3) {
  result <- vector("list", nrow(data))
  keep <- 0L
  for (i in seq_len(nrow(data))) {
    start <- data$midpoint[i]
    end <- start + eclipse_duration_min * 60
    subset <- data %>% filter(midpoint >= start, midpoint <= end)
    if (nrow(subset) >= minimum_dives) {
      keep <- keep + 1L
      result[[keep]] <- data.frame(
        start = start,
        end = end,
        n_dives = nrow(subset),
        mean_residual = mean(subset$adjusted_duration_residual),
        overlaps_eclipse = start <= eclipse_end & end >= eclipse_start
      )
    }
  }
  bind_rows(result[seq_len(keep)])
}

observed_windows <- build_windows(window_analysis, minimum_dives = 3)
control_windows <- observed_windows %>% filter(!overlaps_eclipse)

window_threshold_table <- bind_rows(lapply(window_minimum_dives, function(minimum_dives) {
  controls <- build_windows(window_analysis, minimum_dives) %>% filter(!overlaps_eclipse)
  data.frame(
    comparison = paste0("All controls with at least ", minimum_dives, " dives"),
    eclipse_mean_residual_min = actual_mean_residual,
    n_control_windows = nrow(controls),
    n_equal_or_more_negative = sum(controls$mean_residual <= actual_mean_residual),
    proportion_equal_or_more_negative = mean(controls$mean_residual <= actual_mean_residual),
    median_control_dive_count = median(controls$n_dives),
    control_dive_count_q25 = quantile(controls$n_dives, 0.25),
    control_dive_count_q75 = quantile(controls$n_dives, 0.75)
  )
}))

window_count_rules <- list(
  "Controls with eclipse dive count +/- 3" = function(x) abs(x - eclipse_n_dives) <= 3,
  "Controls with eclipse dive count +/- 2" = function(x) abs(x - eclipse_n_dives) <= 2,
  "Controls with exact eclipse dive count" = function(x) x == eclipse_n_dives
)

window_count_table <- bind_rows(Map(
  function(rule, label) {
    controls <- control_windows[rule(control_windows$n_dives), , drop = FALSE]
    data.frame(
      comparison = label,
      eclipse_mean_residual_min = actual_mean_residual,
      n_control_windows = nrow(controls),
      n_equal_or_more_negative = sum(controls$mean_residual <= actual_mean_residual),
      proportion_equal_or_more_negative = if (nrow(controls) == 0) NA_real_ else
        mean(controls$mean_residual <= actual_mean_residual),
      median_control_dive_count = if (nrow(controls) == 0) NA_real_ else median(controls$n_dives),
      control_dive_count_q25 = if (nrow(controls) == 0) NA_real_ else quantile(controls$n_dives, 0.25),
      control_dive_count_q75 = if (nrow(controls) == 0) NA_real_ else quantile(controls$n_dives, 0.75)
    )
  },
  window_count_rules,
  names(window_count_rules)
))

duration_window_table <- bind_rows(window_threshold_table, window_count_table)

run_moving_block_bootstrap <- function(values, target_n, block_n, B, observed_value) {
  set.seed(20260822 + block_n)
  starts <- seq_len(length(values) - block_n + 1)
  draws <- replicate(B, {
    pseudo <- numeric(0)
    while (length(pseudo) < target_n) {
      start <- sample(starts, 1)
      pseudo <- c(pseudo, values[start:(start + block_n - 1)])
    }
    mean(pseudo[seq_len(target_n)])
  })
  data.frame(
    block_dives = block_n,
    bootstrap_iterations = B,
    n_equal_or_more_negative = sum(draws <= observed_value),
    proportion_equal_or_more_negative = mean(draws <= observed_value),
    lower95 = quantile(draws, 0.025),
    upper95 = quantile(draws, 0.975)
  )
}

if (run_block_bootstrap) {
  residual_series <- window_analysis$adjusted_duration_residual
  block_bootstrap_table <- bind_rows(lapply(bootstrap_block_dives, function(block_n) {
    run_moving_block_bootstrap(residual_series, eclipse_n_dives, block_n,
                               bootstrap_replicates, actual_mean_residual)
  }))
} else {
  block_bootstrap_table <- data.frame(note = "Not run: set run_block_bootstrap = TRUE.")
}

# =============================================================================
# 9. DEPTH-SEQUENCE ANALYSIS
# =============================================================================

sequence_metric <- function(data, flank_n) {
  if (nrow(data) < 2 * flank_n + 1) return(NA_real_)
  first_median <- median(data$dive_depth_model[seq_len(flank_n)])
  last_median <- median(tail(data$dive_depth_model, flank_n))
  reference_depth <- min(first_median, last_median)
  interior <- data$dive_depth_model[(flank_n + 1):(nrow(data) - flank_n)]
  reference_depth - min(interior)
}

build_depth_sequences <- function(data, flank_n) {
  sequences <- vector("list", nrow(data))
  keep <- 0L
  for (i in seq_len(nrow(data))) {
    start <- data$midpoint[i]
    end <- start + eclipse_duration_min * 60
    if (start <= eclipse_end && end >= eclipse_start) next
    subset <- data %>% filter(midpoint >= start, midpoint <= end)
    if (nrow(subset) < 7) next
    value <- sequence_metric(subset, flank_n)
    if (is.finite(value)) {
      keep <- keep + 1L
      sequences[[keep]] <- data.frame(start = start, n_dives = nrow(subset), excursion_m = value)
    }
  }
  bind_rows(sequences[seq_len(keep)])
}

depth_sequence_results <- list()
for (flank_n in 1:3) {
  control_sequences <- build_depth_sequences(window_analysis, flank_n)
  eclipse_excursion <- sequence_metric(actual_eclipse_window, flank_n)
  depth_sequence_results[[as.character(flank_n)]] <- list(
    controls = control_sequences,
    eclipse_excursion = eclipse_excursion
  )
}

depth_sequence_table <- bind_rows(lapply(1:3, function(flank_n) {
  x <- depth_sequence_results[[as.character(flank_n)]]
  data.frame(
    flank_n = flank_n,
    eclipse_excursion_m = x$eclipse_excursion,
    n_control_sequences = nrow(x$controls),
    n_equal_or_larger = sum(x$controls$excursion_m >= x$eclipse_excursion),
    proportion_equal_or_larger = mean(x$controls$excursion_m >= x$eclipse_excursion),
    percentile = 100 * mean(x$controls$excursion_m <= x$eclipse_excursion)
  )
}))

# =============================================================================
# 10. EXPLORATORY OBSCURATION MODELS
# =============================================================================

circle_overlap_fraction <- function(distance) {
  distance <- pmin(pmax(distance, 0), 2)
  (2 * acos(distance / 2) - 0.5 * distance * sqrt(pmax(4 - distance^2, 0))) / pi
}

minimum_disc_separation <- uniroot(
  function(distance) circle_overlap_fraction(distance) - maximum_obscuration,
  interval = c(0, 2)
)$root

# A straight-line relative trajectory is centred halfway between the contacts;
# its maximum differs from the reported maximum time by only 70 s.
geometry_midpoint <- eclipse_start + (eclipse_end - eclipse_start) / 2
geometry_half_duration_min <- as.numeric(difftime(eclipse_end, eclipse_start, units = "mins")) / 2
relative_speed <- sqrt(4 - minimum_disc_separation^2) / geometry_half_duration_min

obscuration_at_time <- function(times) {
  minutes_from_midpoint <- as.numeric(difftime(times, geometry_midpoint, units = "mins"))
  separation <- sqrt(minimum_disc_separation^2 + (relative_speed * minutes_from_midpoint)^2)
  value <- circle_overlap_fraction(separation)
  value[times < eclipse_start | times > eclipse_end] <- 0
  value
}

mean_dive_obscuration <- function(start, end) {
  overlap_start <- max(start, eclipse_start)
  overlap_end <- min(end, eclipse_end)
  if (overlap_end <= overlap_start) return(0)
  times <- seq(overlap_start, overlap_end, by = obscuration_time_step_seconds)
  if (tail(times, 1) != overlap_end) times <- c(times, overlap_end)
  mean(obscuration_at_time(times))
}

obscuration_analysis <- analysis
obscuration_analysis$mean_obscuration <- mapply(
  mean_dive_obscuration, obscuration_analysis$start_time, obscuration_analysis$end_time
)

m_obscuration_duration <- fit_gamm(
  duration_min ~ mean_obscuration + s(solar_elev_deg, k = principal_k, bs = "tp") +
    s(dive_depth_model, k = principal_k, bs = "tp"),
  obscuration_analysis, duration_correlation_type, method = "REML"
)
m_obscuration_depth <- fit_gamm(
  dive_depth_model ~ mean_obscuration + s(solar_elev_deg, k = principal_k, bs = "tp"),
  obscuration_analysis, depth_correlation_type, method = "REML"
)

extract_obscuration_effect <- function(model, response, unit) {
  tab <- summary(model$gam)$p.table
  estimate <- tab["mean_obscuration", "Estimate"]
  se <- tab["mean_obscuration", "Std. Error"]
  data.frame(
    response = response, unit = unit,
    estimate_per_full_obscuration = estimate,
    lower95 = estimate - 1.96 * se,
    upper95 = estimate + 1.96 * se,
    p = tab["mean_obscuration", ncol(tab)]
  )
}

obscuration_table <- bind_rows(
  extract_obscuration_effect(m_obscuration_duration, "Dive duration", "min"),
  extract_obscuration_effect(m_obscuration_depth, "Dive depth", "m")
)

# =============================================================================
# 11. TABLES, AUDIT DATA, MODELS, AND REPORT
# =============================================================================

table_s1 <- data.frame(
  quantity = c(
    "Raw dive records from deployment onwards",
    "Dive records with daily high-quality Argos position",
    "Dive records from complete recording days",
    "Eclipse dives defined by midpoint",
    "Eclipse dives defined by any temporal overlap",
    "Median duration-bin width (s)",
    "Median depth-bin width (m)"
  ),
  value = c(
    nrow(raw_dive_records), nrow(dives),
    sum(dives$date != deployment_date & dives$date != last_incomplete_date),
    sum(dives$eclipse_midpoint), sum(dives$eclipse_overlap),
    median(raw_dive_records$duration_bin_width_sec, na.rm = TRUE),
    median(raw_dive_records$depth_bin_width_m, na.rm = TRUE)
  )
)

table_s2 <- bind_rows(principal_duration_effect, principal_depth_effect) %>%
  transmute(
    analysis = paste("Principal", response, "model"),
    unit, n_dives = nrow(analysis), estimate, lower95, upper95, p
  )

write_table(table_s1, "Table_S1_data_processing_summary.csv")
write_table(table_s2, "Table_S2_principal_model_results.csv")
write_table(temporal_table, "Table_S3_temporal_correlation_selection.csv")
write_table(duration_sensitivity_table, "Table_S4_duration_model_sensitivities.csv")
write_table(counterfactual_table, "Table_S5_counterfactual_duration.csv")
write_table(duration_window_table, "Table_S6_duration_window_analysis.csv")
write_table(depth_sequence_table, "Table_S7_depth_sequence_analysis.csv")
write_table(obscuration_table, "Table_S8_obscuration_analysis.csv")

write_table(duration_consistency, "duration_consistency.csv", audit_dir)
write_table(depth_resolution, "depth_resolution.csv", audit_dir)
write_table(eclipse_definition_summary, "eclipse_definition_summary.csv", audit_dir)
write_table(common_support_summary, "principal_common_support.csv", audit_dir)
write_table(local_support_table, "principal_local_control_density.csv", audit_dir)
write_table(nearest_control_table, "principal_nearest_controls.csv", audit_dir)
write_table(block_bootstrap_table, "moving_block_bootstrap.csv", audit_dir)

write.csv(dives, file.path(data_dir, "complete_processed_dive_record.csv"), row.names = FALSE)
write.csv(analysis, file.path(data_dir, "principal_analysis_dataset.csv"), row.names = FALSE)
write.csv(window_analysis, file.path(data_dir, "complete_day_window_dataset.csv"), row.names = FALSE)
write.csv(obscuration_analysis, file.path(data_dir, "obscuration_analysis_dataset.csv"), row.names = FALSE)

saveRDS(
  list(
    baseline_duration = m_base_duration,
    baseline_depth = m_base_depth,
    principal_duration = m_duration,
    principal_depth = m_depth,
    window_reference_duration = m_window_reference,
    obscuration_duration = m_obscuration_duration,
    obscuration_depth = m_obscuration_depth,
    duration_correlation_type = duration_correlation_type,
    depth_correlation_type = depth_correlation_type
  ),
  file.path(model_dir, "supplementary_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

report_file <- file.path(output_dir, "FINAL_ANALYSIS_REPORT.txt")
report_lines <- capture.output({
  cat("HUMPBACK WHALE SOLAR ECLIPSE — FINAL ANALYSIS\n")
  cat("================================================\n\n")
  cat("Temporal structure selected by ML AIC:\n")
  cat("  Dive duration:", temporal_labels[duration_correlation_type], "\n")
  cat("  Dive depth:", temporal_labels[depth_correlation_type], "\n\n")
  cat("Principal duration model (REML):\n")
  print(principal_duration_effect, row.names = FALSE)
  cat("\nPrincipal depth model (REML):\n")
  print(principal_depth_effect, row.names = FALSE)
  cat("\nTemporal-correlation comparison (ML):\n")
  print(temporal_table, row.names = FALSE)
  cat("\nCounterfactual duration effect:\n")
  print(counterfactual_table, row.names = FALSE)
  cat("\nDuration-window comparison:\n")
  print(duration_window_table, row.names = FALSE)
  cat("\nDepth-sequence comparison:\n")
  print(depth_sequence_table, row.names = FALSE)
  cat("\nObscuration models:\n")
  print(obscuration_table, row.names = FALSE)
})
writeLines(report_lines, report_file)

# =============================================================================
# 12. MANUSCRIPT FIGURE 4
# =============================================================================

# This recreates the current manuscript layout directly from the objects fitted
# above. It is deliberately produced here, rather than by a separate script, so
# every annotation and line is tied to the final models and output tables.

font_main <- "Helvetica"
col_black <- "#171717"
col_axis <- "#3F3F3F"
col_points <- "#AFAFAF"
col_points_dark <- "#747474"
col_band <- "#D9D9D9"
col_hist <- "#BCBCBC"
col_eclipse <- "#E76F00"
col_eclipse_fill <- "#F4D9CA"

theme_figure4 <- theme_classic(base_family = font_main, base_size = 10.5) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = col_axis, linewidth = 0.45),
    axis.ticks = element_line(colour = col_axis, linewidth = 0.35),
    axis.ticks.length = unit(2.5, "pt"),
    axis.text = element_text(family = font_main, colour = col_axis, size = 9),
    axis.title = element_text(family = font_main, colour = col_black, size = 10.2),
    plot.margin = margin(5, 8, 5, 7),
    legend.position = "none"
  )

panel_title_grob <- function(letter, title) {
  grobTree(
    rectGrob(
      x = unit(0, "npc"), y = unit(0.52, "npc"),
      width = unit(0.038, "npc"), height = unit(0.54, "npc"),
      just = c("left", "center"), gp = gpar(fill = col_black, col = NA)
    ),
    textGrob(
      letter, x = unit(0.019, "npc"), y = unit(0.52, "npc"),
      gp = gpar(fontfamily = font_main, fontsize = 10.5, fontface = "bold", col = "white")
    ),
    textGrob(
      title, x = unit(0.054, "npc"), y = unit(0.52, "npc"),
      just = c("left", "center"),
      gp = gpar(fontfamily = font_main, fontsize = 11.4, col = "#111111")
    )
  )
}

with_panel_title <- function(plot, letter, title, title_height = 0.09) {
  arrangeGrob(
    panel_title_grob(letter, title), ggplotGrob(plot),
    ncol = 1, heights = c(title_height, 1)
  )
}

plot_dives <- dives %>%
  filter(is.finite(solar_elev_deg), is.finite(duration_min), is.finite(dive_depth_m)) %>%
  arrange(midpoint)

solar_time_grid <- data.frame(
  datetime = seq(min(plot_dives$midpoint), max(plot_dives$midpoint), by = "2 min")
) %>%
  mutate(date = as.Date(datetime, tz = "UTC")) %>%
  left_join(daily_position, by = "date")

if (any(!is.finite(solar_time_grid$daily_lat) | !is.finite(solar_time_grid$daily_lon))) {
  stop("At least one date in the Figure 4 solar grid lacks a daily Argos position.")
}

solar_grid_position <- suncalc::getSunlightPosition(data = data.frame(
  date = solar_time_grid$datetime,
  lat = solar_time_grid$daily_lat,
  lon = solar_time_grid$daily_lon
))
solar_time_grid$solar_elev_deg <- solar_grid_position$altitude * 180 / pi

p_solar <- ggplot(solar_time_grid, aes(datetime, solar_elev_deg)) +
  annotate("rect", xmin = eclipse_start, xmax = eclipse_end, ymin = -Inf, ymax = Inf,
           fill = col_eclipse_fill, alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#A0A0A0", linewidth = 0.35) +
  geom_area(fill = "#EEEEEE", alpha = 0.75) +
  geom_line(colour = col_black, linewidth = 0.68, lineend = "round") +
  geom_vline(xintercept = eclipse_max, colour = col_eclipse, linewidth = 0.85) +
  annotate("text", x = eclipse_max, y = Inf, label = "Solar eclipse",
           colour = col_eclipse, family = font_main, size = 3, hjust = -0.1, vjust = 1.2) +
  scale_x_datetime(date_breaks = "2 days", date_labels = "%d Aug",
                   expand = expansion(mult = c(0.004, 0.004))) +
  scale_y_continuous(breaks = c(-10, 0, 20, 40),
                     expand = expansion(mult = c(0.04, 0.09))) +
  labs(title = "Solar elevation throughout the behavioural record", x = NULL,
       y = "Solar elevation (°)") +
  theme_figure4 +
  theme(
    plot.title = element_text(family = font_main, size = 10.6, margin = margin(b = 3)),
    axis.title.y = element_text(size = 9.2), axis.text = element_text(size = 8.4),
    plot.margin = margin(2, 8, 1, 7)
  )

prediction_grid <- function(model, x_name, x_values) {
  new_data <- setNames(data.frame(x_values), x_name)
  prediction <- predict(model$gam, newdata = new_data, type = "link", se.fit = TRUE)
  new_data$fit <- as.numeric(prediction$fit)
  new_data$lower <- new_data$fit - 1.96 * as.numeric(prediction$se.fit)
  new_data$upper <- new_data$fit + 1.96 * as.numeric(prediction$se.fit)
  new_data
}

solar_grid_A <- prediction_grid(
  m_base_duration, "solar_elev_deg",
  seq(min(plot_dives$solar_elev_deg), max(plot_dives$solar_elev_deg), length.out = 500)
)
solar_grid_B <- prediction_grid(
  m_base_depth, "solar_elev_deg",
  seq(min(plot_dives$solar_elev_deg), max(plot_dives$solar_elev_deg), length.out = 500)
)

duration_solar_p <- summary(m_base_duration$gam)$s.table[1, "p-value"]
depth_solar_p <- summary(m_base_depth$gam)$s.table[1, "p-value"]

p_A <- ggplot() +
  annotate("rect", xmin = eclipse_solar_min, xmax = eclipse_solar_max, ymin = -Inf, ymax = Inf,
           fill = col_eclipse_fill, alpha = 0.48) +
  geom_ribbon(data = solar_grid_A, aes(solar_elev_deg, ymin = lower, ymax = upper),
              fill = col_band, alpha = 0.65) +
  geom_point(data = filter(plot_dives, !eclipse), aes(solar_elev_deg, duration_min),
             colour = col_points, alpha = 0.35, size = 1.05) +
  geom_line(data = solar_grid_A, aes(solar_elev_deg, fit), colour = col_black, linewidth = 0.9) +
  geom_point(data = filter(plot_dives, eclipse), aes(solar_elev_deg, duration_min),
             colour = col_eclipse, size = 1.85) +
  geom_vline(xintercept = c(eclipse_solar_min, eclipse_solar_max), linetype = "dashed",
             colour = col_eclipse, linewidth = 0.35) +
  annotate("text", x = mean(c(eclipse_solar_min, eclipse_solar_max)), y = Inf,
           label = "Eclipse solar range", family = font_main, colour = col_eclipse,
           size = 2.75, hjust = 0.5, vjust = 1.3) +
  annotate("text", x = -Inf, y = Inf, label = format_p(duration_solar_p),
           family = font_main, colour = "#555555", size = 2.9, hjust = -0.08, vjust = 1.3) +
  labs(x = "Solar elevation (°)", y = "Dive duration (min)") +
  theme_figure4

p_B <- ggplot() +
  annotate("rect", xmin = eclipse_solar_min, xmax = eclipse_solar_max, ymin = -Inf, ymax = Inf,
           fill = col_eclipse_fill, alpha = 0.48) +
  geom_ribbon(data = solar_grid_B, aes(solar_elev_deg, ymin = lower, ymax = upper),
              fill = col_band, alpha = 0.65) +
  geom_point(data = filter(plot_dives, !eclipse), aes(solar_elev_deg, dive_depth_m),
             colour = col_points, alpha = 0.35, size = 1.05) +
  geom_line(data = solar_grid_B, aes(solar_elev_deg, fit), colour = col_black, linewidth = 0.9) +
  geom_point(data = filter(plot_dives, eclipse), aes(solar_elev_deg, dive_depth_m),
             colour = col_eclipse, size = 1.85) +
  geom_vline(xintercept = c(eclipse_solar_min, eclipse_solar_max), linetype = "dashed",
             colour = col_eclipse, linewidth = 0.35) +
  annotate("text", x = mean(c(eclipse_solar_min, eclipse_solar_max)), y = 0,
           label = "Eclipse solar range", family = font_main, colour = col_eclipse,
           size = 2.75, hjust = 0.5, vjust = -0.75) +
  annotate("text", x = -Inf, y = 0, label = format_p(depth_solar_p),
           family = font_main, colour = "#555555", size = 2.9, hjust = -0.08, vjust = -1.3) +
  scale_y_reverse(breaks = c(0, 100, 200, 300, 400)) +
  labs(x = "Solar elevation (°)", y = "Dive depth (m)") +
  theme_figure4

# Panel C: conditional duration predictions across dive depth.
depth_grid_C <- seq(min(analysis$dive_depth_model), max(analysis$dive_depth_model), length.out = 400)
solar_reference_C <- mean(actual_eclipse$solar_elev_deg)
prediction_control_C <- data.frame(
  solar_elev_deg = solar_reference_C,
  dive_depth_model = depth_grid_C,
  eclipse_f = factor("control", levels = c("control", "eclipse"))
)
prediction_eclipse_C <- prediction_control_C
prediction_eclipse_C$eclipse_f <- factor("eclipse", levels = c("control", "eclipse"))

add_prediction <- function(data) {
  prediction <- predict(m_duration$gam, newdata = data, type = "link", se.fit = TRUE)
  data$fit <- as.numeric(prediction$fit)
  data$lower <- pmax(0, data$fit - 1.96 * as.numeric(prediction$se.fit))
  data$upper <- data$fit + 1.96 * as.numeric(prediction$se.fit)
  data
}
prediction_control_C <- add_prediction(prediction_control_C)
prediction_eclipse_C <- add_prediction(prediction_eclipse_C)

p_C_top <- ggplot() +
  geom_density(data = filter(analysis, !eclipse_model), aes(dive_depth_model),
               fill = "#BFBFBF", colour = NA, alpha = 0.68) +
  geom_density(data = actual_eclipse, aes(dive_depth_model),
               fill = col_eclipse, colour = NA, alpha = 0.52, adjust = 0.85) +
  scale_x_continuous(limits = range(analysis$dive_depth_model), expand = c(0, 0)) +
  theme_void(base_family = font_main) + theme(plot.margin = margin(0, 0, -2, 0))

p_C_right <- ggplot() +
  geom_density(data = filter(analysis, !eclipse_model), aes(duration_min),
               fill = "#BFBFBF", colour = NA, alpha = 0.68) +
  geom_density(data = actual_eclipse, aes(duration_min),
               fill = col_eclipse, colour = NA, alpha = 0.52, adjust = 0.85) +
  scale_x_continuous(limits = c(0, max(analysis$duration_min)), expand = c(0, 0)) +
  coord_flip() + theme_void(base_family = font_main) +
  theme(plot.margin = margin(0, 0, 0, -2))

effect_label <- paste0(
  "Eclipse effect = ", sprintf("%.2f", principal_duration_effect$estimate), " min\n",
  "95% CI: ", sprintf("%.2f", principal_duration_effect$lower95), " to ",
  sprintf("%.2f", principal_duration_effect$upper95), "\n",
  format_p(principal_duration_effect$p)
)

p_C_main <- ggplot() +
  geom_ribbon(data = prediction_control_C,
              aes(dive_depth_model, ymin = lower, ymax = upper), fill = col_band, alpha = 0.62) +
  geom_ribbon(data = prediction_eclipse_C,
              aes(dive_depth_model, ymin = lower, ymax = upper), fill = col_eclipse_fill, alpha = 0.72) +
  geom_point(data = filter(analysis, !eclipse_model), aes(dive_depth_model, duration_min),
             colour = col_points, alpha = 0.24, size = 0.95) +
  geom_line(data = prediction_control_C, aes(dive_depth_model, fit),
            colour = col_black, linewidth = 0.9) +
  geom_line(data = prediction_eclipse_C, aes(dive_depth_model, fit),
            colour = col_eclipse, linewidth = 0.9) +
  geom_point(data = actual_eclipse, aes(dive_depth_model, duration_min),
             colour = col_eclipse, size = 1.75) +
  annotate("text", x = Inf, y = Inf, label = effect_label, family = font_main,
           colour = col_black, size = 2.85, lineheight = 1.08, hjust = 1.03, vjust = 1.12) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Dive depth (m)", y = "Dive duration (min)") +
  theme_figure4 + theme(plot.margin = margin(0, 0, 4, 7))

g_C_top <- ggplotGrob(p_C_top)
g_C_main <- ggplotGrob(p_C_main)
g_C_right <- ggplotGrob(p_C_right)
max_C_widths <- unit.pmax(g_C_top$widths, g_C_main$widths)
g_C_top$widths <- max_C_widths
g_C_main$widths <- max_C_widths
g_C_body <- arrangeGrob(g_C_top, nullGrob(), g_C_main, g_C_right,
                        nrow = 2, ncol = 2, widths = c(1, 0.15), heights = c(0.20, 1))
g_C <- arrangeGrob(
  panel_title_grob("C", "Predicted duration under non-eclipse and eclipse conditions"),
  g_C_body, ncol = 1, heights = c(0.09, 1)
)

# Panel D: empirical 125-min contextual comparison.
n_windows_D <- nrow(control_windows)
n_extreme_D <- sum(control_windows$mean_residual <= actual_mean_residual)
pct_extreme_D <- 100 * n_extreme_D / n_windows_D
p_D <- ggplot(control_windows, aes(mean_residual)) +
  annotate("rect", xmin = -Inf, xmax = actual_mean_residual, ymin = -Inf, ymax = Inf,
           fill = col_eclipse_fill, alpha = 0.32) +
  geom_histogram(binwidth = 0.25, boundary = 0, fill = col_hist, colour = "white", linewidth = 0.32) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#555555", linewidth = 0.58) +
  geom_vline(xintercept = actual_mean_residual, colour = col_eclipse, linewidth = 1.15) +
  annotate("text", x = actual_mean_residual, y = Inf,
           label = paste0("Eclipse = ", sprintf("%.2f", actual_mean_residual), " min"),
           colour = col_eclipse, family = font_main, size = 3, hjust = -0.08, vjust = 1.25) +
  annotate("text", x = Inf, y = Inf,
           label = paste0(n_extreme_D, " / ", format(n_windows_D, big.mark = ","), " windows\n",
                          sprintf("%.2f", pct_extreme_D), "% equally or more negative"),
           family = font_main, size = 3, lineheight = 1.08, hjust = 1.03, vjust = 1.18) +
  labs(x = "Mean observed − expected dive duration\nper 125-min window (min)",
       y = "Number of windows") +
  theme_figure4

# Global legend and final Figure 4 grid.
p_legend <- ggplot() +
  annotate("point", x = 0.05, y = 0.5, colour = col_points_dark, size = 2.3) +
  annotate("text", x = 0.075, y = 0.5, label = "Non-eclipse dives", family = font_main, size = 3.05, hjust = 0) +
  annotate("point", x = 0.265, y = 0.5, colour = col_eclipse, size = 2.5) +
  annotate("text", x = 0.290, y = 0.5, label = "Eclipse dives", family = font_main, size = 3.05, hjust = 0) +
  annotate("segment", x = 0.455, xend = 0.500, y = 0.5, yend = 0.5, colour = col_black, linewidth = 0.9) +
  annotate("text", x = 0.515, y = 0.5, label = "Expected / non-eclipse prediction", family = font_main, size = 3.05, hjust = 0) +
  annotate("segment", x = 0.785, xend = 0.830, y = 0.5, yend = 0.5, colour = col_eclipse, linewidth = 0.9) +
  annotate("text", x = 0.845, y = 0.5, label = "Eclipse prediction", family = font_main, size = 3.05, hjust = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

g_A <- with_panel_title(p_A, "A", "Solar-associated variation in dive duration")
g_B <- with_panel_title(p_B, "B", "Solar-associated variation in dive depth")
g_D <- with_panel_title(p_D, "D", "Eclipse relative to equivalent periods in the behavioural record")
g_solar <- ggplotGrob(p_solar)
g_legend <- ggplotGrob(p_legend)

draw_figure4 <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 4, ncol = 2,
    heights = unit(c(0.63, 0.17, 1.30, 1.90), "null"),
    widths = unit(c(1.07, 0.93), "null")
  )
  pushViewport(viewport(layout = layout))
  positions <- list(
    list(grob = g_solar, row = 1, col = 1:2),
    list(grob = g_legend, row = 2, col = 1:2),
    list(grob = g_A, row = 3, col = 1),
    list(grob = g_B, row = 3, col = 2),
    list(grob = g_C, row = 4, col = 1),
    list(grob = g_D, row = 4, col = 2)
  )
  for (position in positions) {
    pushViewport(viewport(layout.pos.row = position$row, layout.pos.col = position$col))
    grid.draw(position$grob)
    popViewport()
  }
  popViewport()
}

figure4_svg <- file.path(figure_dir, "Figure4_final.svg")
svglite::svglite(figure4_svg, width = 13, height = 9.8, bg = "white")
draw_figure4()
dev.off()

figure4_pdf <- file.path(figure_dir, "Figure4_final.pdf")
cairo_pdf(figure4_pdf, width = 13, height = 9.8, bg = "white")
draw_figure4()
dev.off()

figure4_png <- file.path(figure_dir, "Figure4_final.png")
png(figure4_png, width = 13, height = 9.8, units = "in", res = 600, bg = "white")
draw_figure4()
dev.off()

# =============================================================================
# 13. SUPPLEMENTARY FIGURES
# =============================================================================

duration_residuals <- residuals(m_duration$lme, type = "normalized")
duration_fitted <- fitted(m_duration$lme)
qq <- qqnorm(duration_residuals, plot.it = FALSE)
acf_values <- acf(duration_residuals, lag.max = 30, plot = FALSE)
acf_data <- data.frame(lag = acf_values$lag[-1], acf = as.numeric(acf_values$acf[-1]))
acf_bound <- 1.96 / sqrt(length(duration_residuals))

p_residual_fitted <- ggplot(data.frame(fitted = duration_fitted, residual = duration_residuals),
                            aes(fitted, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, size = 1.2) +
  labs(x = "Fitted duration", y = "Normalised residual") + theme_classic()

p_qq <- ggplot(data.frame(theoretical = qq$x, sample = qq$y), aes(theoretical, sample)) +
  geom_point(alpha = 0.65, size = 1.2) + geom_abline(linetype = "dashed", colour = "grey50") +
  labs(x = "Theoretical quantiles", y = "Normalised residual quantiles") + theme_classic()

p_acf <- ggplot(acf_data, aes(lag, acf)) +
  geom_hline(yintercept = c(-acf_bound, acf_bound), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0) + geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  labs(x = "Lag (dives)", y = "ACF") + theme_classic()

# Breaks stop exactly at the selected support limits, avoiding the visual
# overhang of the terminal histogram bin beyond the orange selection line.
p_support <- ggplot(support_control, aes(solar_elev_deg)) +
  geom_histogram(
    breaks = seq(eclipse_solar_min, eclipse_solar_max, length.out = 21),
    fill = "grey80", colour = "white"
  ) +
  geom_vline(xintercept = c(eclipse_solar_min, eclipse_solar_max),
             colour = col_eclipse, linewidth = 0.8) +
  labs(x = "Solar elevation (°)", y = "Control dives") + theme_classic()

figure_s1 <- arrangeGrob(
  p_residual_fitted, p_qq, p_acf, p_support,
  ncol = 2, top = textGrob("Principal duration-model diagnostics", gp = gpar(fontsize = 14))
)
save_figure(figure_s1, "Figure_S1_principal_model_diagnostics", width = 8, height = 6.5)

primary_depth_sequences <- depth_sequence_results[["2"]]$controls
primary_eclipse_excursion <- depth_sequence_results[["2"]]$eclipse_excursion
p_depth_context <- ggplot(window_analysis, aes(midpoint, dive_depth_model)) +
  annotate("rect", xmin = eclipse_start, xmax = eclipse_end, ymin = -Inf, ymax = Inf,
           fill = col_eclipse_fill, alpha = 0.6) +
  geom_line(colour = "grey60", linewidth = 0.35) +
  geom_point(aes(colour = eclipse_model), size = 1.15) +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = col_eclipse), guide = "none") +
  scale_y_reverse() +
  labs(x = "Date (UTC)", y = "Dive depth (m)") + theme_classic()

p_excursion <- ggplot(primary_depth_sequences, aes(excursion_m)) +
  geom_histogram(binwidth = 10, fill = "grey75", colour = "white") +
  geom_vline(xintercept = primary_eclipse_excursion, colour = col_eclipse, linewidth = 1) +
  labs(x = "Shallow excursion (m)", y = "Control sequences") + theme_classic()

figure_s2 <- arrangeGrob(
  p_depth_context, p_excursion, ncol = 1,
  top = textGrob("Depth-sequence context", gp = gpar(fontsize = 14))
)
save_figure(figure_s2, "Figure_S2_depth_sequence_context", width = 9, height = 7)

# =============================================================================
# 14. FINAL CONSOLE SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("FINAL ECLIPSE ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Outputs:", output_dir, "\n")
cat("Duration correlation:", temporal_labels[duration_correlation_type], "\n")
cat("Depth correlation:", temporal_labels[depth_correlation_type], "\n")
cat("Figure 4:", figure4_svg, "\n")
cat("Tables S1–S8:", table_dir, "\n")
