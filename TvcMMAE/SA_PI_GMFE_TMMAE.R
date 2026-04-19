library(ospsuite)
library(ggplot2)

rm(list = ls())

options(width = 140)
options(ospsuite.plots.watermarkEnabled = FALSE)

# Paths -------------------------------------------------------------------

sim_path <- "~/Desktop/Thesis/TvcMMAE/Mouse TvcMMAE i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TvcMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/TvcMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# Plot options ------------------------------------------------------------

include_observed <- TRUE
plot_title <- "Plasma TvcMMAE Concentration-Time Profiles"
x_axis_label <- "Time [h]"
y_axis_label <- "Concentration [\u03bcmol/l]"
legend_position <- "right"
y_scale <- "log"

obs_label <- "Observed data"
obs_color <- "orange"
obs_shape <- 16
obs_size <- 2.2

# Observed data -----------------------------------------------------------

obs_data <- loadDataSetFromPKML(filePath = dataset_path)
print(obs_data)

# Helper functions --------------------------------------------------------

to_df <- function(x, pk_name) {
  if (length(x) == 0) {
    return(
      data.frame(
        pk_parameter = character(0),
        parameter = character(0),
        sensitivity = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    pk_parameter = rep(pk_name, length(x)),
    parameter = vapply(x, \(y) y$parameterName, character(1)),
    sensitivity = vapply(x, \(y) y$value, numeric(1)),
    stringsAsFactors = FALSE
  )
}

apply_changes <- function(
    sim,
    fcrn_endo_kd = NULL,
    tmmae_her2_kd = NULL,
    herref = NULL,
    dose = NULL,
    fu = NULL
) {
  if (!is.null(fu)) {
    fu_plasma <- getParameter(
      "TvcMMAE|Fraction unbound (plasma, reference value)",
      sim
    )
    setParameterValues(fu_plasma, fu)
  }
  
  if (!is.null(tmmae_her2_kd)) {
    her2_kd <- getParameter("TvcMMAE-HER2-Chang et al.|Kd", sim)
    setParameterValues(her2_kd, tmmae_her2_kd)
  }
  
  if (!is.null(fcrn_endo_kd)) {
    fcrn_kd <- getParameter(
      "**|TvcMMAE|Kd (FcRn) in Endosomal Space",
      sim
    )
    setParameterValues(fcrn_kd, fcrn_endo_kd)
  }
  
  if (!is.null(dose)) {
    dose_par <- getParameter(
      "Events|10 mg/kg|Application_1|ProtocolSchemaItem|DosePerBodyWeight",
      sim
    )
    setParameterValues(dose_par, dose)
  }
  
  if (!is.null(herref)) {
    herref_par <- getParameter("HER2|Reference concentration", sim)
    setParameterValues(herref_par, herref)
  }
  
  sim$solver$absTol <- 1e-12
  sim$solver$relTol <- 1e-8
  
  invisible(sim)
}

get_last_obs_time_h <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  max(obs_df$xValues, na.rm = TRUE)
}

get_obs_auc_to_last_obs <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  obs_df <- obs_df[order(obs_df$xValues), ]
  
  x <- obs_df$xValues
  y <- obs_df$yValues
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  auc_h <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  auc_min <- auc_h * 60
  
  auc_min
}

define_auc_to_last_obs_parameter <- function(end_time_h) {
  removeAllUserDefinedPKParameters()
  
  auc_to_last_obs <- addUserDefinedPKParameter(
    name = "AUC_to_lastObs",
    standardPKParameter = StandardPKParameter$AUC_tEnd
  )
  
  auc_to_last_obs$startTime <- 0
  auc_to_last_obs$endTime <- end_time_h * 60
  
  invisible(auc_to_last_obs)
}

get_sim_pk_value <- function(sim_results, output_path, parameter_name) {
  pk <- calculatePKAnalyses(results = sim_results)
  pk_df <- pkAnalysesToDataFrame(pk)
  
  row <- pk_df[
    pk_df$QuantityPath == output_path &
      pk_df$Parameter == parameter_name,
  ]
  
  if (nrow(row) == 0) {
    stop(paste("PK parameter not found:", parameter_name))
  }
  
  row$Value[[1]]
}

fold_error <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  out <- rep(NA_real_, length(pred))
  out[ok] <- 10^(abs(log10(pred[ok] / obs[ok])))
  out
}

build_legend_spec <- function(
    curve_spec,
    include_observed = FALSE,
    obs_label = "Observed data",
    obs_color = "orange",
    obs_shape = 16
) {
  color_values <- stats::setNames(curve_spec$color, curve_spec$label)
  linetype_values <- stats::setNames(curve_spec$linetype, curve_spec$label)
  linewidth_values <- stats::setNames(curve_spec$linewidth, curve_spec$label)
  shape_values <- stats::setNames(rep(NA, nrow(curve_spec)), curve_spec$label)
  
  if (include_observed) {
    color_values <- c(
      color_values,
      stats::setNames(obs_color, obs_label)
    )
    linetype_values <- c(
      linetype_values,
      stats::setNames("blank", obs_label)
    )
    linewidth_values <- c(
      linewidth_values,
      stats::setNames(0, obs_label)
    )
    shape_values <- c(
      shape_values,
      stats::setNames(obs_shape, obs_label)
    )
    legend_order <- c(curve_spec$label, obs_label)
  } else {
    legend_order <- curve_spec$label
  }
  
  list(
    order = legend_order,
    color_values = color_values,
    linetype_values = linetype_values,
    linewidth_values = linewidth_values,
    shape_values = shape_values
  )
}

make_scenario_plot <- function(
    plot_data,
    legend_spec,
    title,
    include_observed = FALSE,
    x_unit = ospUnits$Time$h,
    y_scale = "log",
    x_axis_label = "Time [h]",
    y_axis_label = "Concentration [\u03bcmol/l]",
    legend_position = "right"
) {
  base_plot <- plotTimeProfile(
    plotData = plot_data,
    xUnit = x_unit,
    yScale = y_scale,
    mapping = aes(
      color = name,
      linetype = name,
      linewidth = name
    ),
    observedMapping = if (include_observed) {
      aes(color = name, shape = name)
    } else {
      NULL
    },
    showLegendPerDataset = "all"
  )
  
  base_plot +
    scale_color_manual(
      values = legend_spec$color_values,
      breaks = legend_spec$order
    ) +
    scale_linetype_manual(
      values = legend_spec$linetype_values,
      breaks = legend_spec$order
    ) +
    scale_linewidth_manual(
      values = legend_spec$linewidth_values,
      breaks = legend_spec$order
    ) +
    scale_shape_manual(
      values = legend_spec$shape_values,
      breaks = legend_spec$order
    ) +
    labs(
      title = title,
      x = x_axis_label,
      y = y_axis_label,
      color = NULL,
      linetype = NULL,
      linewidth = NULL,
      shape = NULL
    ) +
    guides(
      linetype = "none",
      linewidth = "none",
      shape = "none",
      color = guide_legend(
        override.aes = list(
          linetype = unname(legend_spec$linetype_values[legend_spec$order]),
          linewidth = unname(legend_spec$linewidth_values[legend_spec$order]),
          shape = unname(legend_spec$shape_values[legend_spec$order])
        )
      )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = legend_position,
      legend.title = element_blank()
    )
}

# Baseline sensitivity analysis -------------------------------------------

sim_sa <- loadSimulation(sim_path)
setOutputs(out_path, sim_sa)

sa_parameter_paths <- c(
  "**|TvcMMAE|Kd (FcRn) in Endosomal Space",
  "TvcMMAE-HER2-Chang et al|Kd",
  "HER2|Reference concentration",
  "TvcMMAE|Fraction unbound (plasma, reference value)"
)

sa <- SensitivityAnalysis$new(
  simulation = sim_sa,
  parameterPaths = sa_parameter_paths,
  variationRange = 0.1,
  numberOfSteps = 2
)

sa_result <- runSensitivityAnalysis(sa)

sens_auc_tend <- sa_result$allPKParameterSensitivitiesFor(
  pkParameterName = "AUC_tEnd",
  outputPath = out_path,
  totalSensitivityThreshold = 1
)

sens_cl <- sa_result$allPKParameterSensitivitiesFor(
  pkParameterName = "CL",
  outputPath = out_path,
  totalSensitivityThreshold = 1
)

sa_df <- rbind(
  to_df(sens_auc_tend, "AUC_tEnd"),
  to_df(sens_cl, "CL")
)

sa_df <- sa_df[order(sa_df$pk_parameter, -abs(sa_df$sensitivity)), ]

cat("\nSensitivity summary:\n")
print(sa_df, row.names = FALSE)

# Scenarios ---------------------------------------------------------------

scenarios <- list(
  list(
    id = "S1",
    label = "FcRnKd: 1.07 \u03bcmol/l, HER2Kd: 1 \u03bcmol/l, HER2 ref: 0.7 \u03bcmol/l",
    fcrn_endo_kd = 1.07,
    tmmae_her2_kd = 1,
    herref = 0.7,
    dose = NULL,
    fu = NULL,
    color = "red",
    linewidth = 0.9,
    linetype = "solid"
  ),
  list(
    id = "S2",
    label = "FcRnKd: 1.07 \u03bcmol/l, HER2Kd: 8 \u03bcmol/l, HER2 ref: 0.25 \u03bcmol/l, fu: 0.5",
    fcrn_endo_kd = 1.07,
    tmmae_her2_kd = 8,
    herref = 0.25,
    dose = NULL,
    fu = 0.5,
    color = "blue",
    linewidth = 0.9,
    linetype = "dashed"
  ),
  list(
    id = "S3",
    label = "FcRnKd: 1.07 \u03bcmol/l, HER2Kd: 10 \u03bcmol/l, HER2 ref: 0.1 \u03bcmol/l, fu: 0.75",
    fcrn_endo_kd = 1.07,
    tmmae_her2_kd = 10,
    herref = 0.1,
    dose = NULL,
    fu = 0.75,
    color = "darkgreen",
    linewidth = 0.9,
    linetype = "dotdash"
  )
)

# Observed reference metric -----------------------------------------------

last_obs_time_h <- get_last_obs_time_h(obs_data)
obs_auc_to_last_obs <- get_obs_auc_to_last_obs(obs_data)

define_auc_to_last_obs_parameter(last_obs_time_h)

cat("\nObserved reference:\n")
cat("Last observed time (h):", round(last_obs_time_h, 3), "\n")
cat("Observed AUC_to_lastObs (\u03bcmol*min/l):", round(obs_auc_to_last_obs, 3), "\n")

# Combined data for plot + scenario PK table ------------------------------

dc <- DataCombined$new()

if (include_observed) {
  obs_data$name <- obs_label
  dc$addDataSets(obs_data)
}

scenario_pk_list <- vector("list", length(scenarios))
scenario_lookup_list <- vector("list", length(scenarios))

for (i in seq_along(scenarios)) {
  scn <- scenarios[[i]]
  
  sim <- loadSimulation(sim_path)
  setOutputs(out_path, sim)
  
  apply_changes(
    sim,
    fcrn_endo_kd = scn$fcrn_endo_kd,
    tmmae_her2_kd = scn$tmmae_her2_kd,
    herref = scn$herref,
    dose = scn$dose,
    fu = scn$fu
  )
  
  res <- runSimulations(sim)[[1]]
  
  sim_auc_to_last_obs <- get_sim_pk_value(res, out_path, "AUC_to_lastObs")
  
  scenario_pk_list[[i]] <- data.frame(
    scenario = scn$id,
    auc_to_last_obs = sim_auc_to_last_obs,
    stringsAsFactors = FALSE
  )
  
  scenario_lookup_list[[i]] <- data.frame(
    scenario = scn$id,
    description = scn$label,
    stringsAsFactors = FALSE
  )
  
  dc$addSimulationResults(
    simulationResults = res,
    quantitiesOrPaths = out_path,
    names = scn$id
  )
}

scenario_pk <- do.call(rbind, scenario_pk_list)
scenario_lookup <- do.call(rbind, scenario_lookup_list)

# Error analysis ----------------------------------------------------------

error_summary <- data.frame(
  scenario = scenario_pk$scenario,
  observed_auc_to_last_obs = round(
    rep(obs_auc_to_last_obs, nrow(scenario_pk)),
    3
  ),
  predicted_auc_to_last_obs = round(scenario_pk$auc_to_last_obs, 3),
  fold_error = round(
    fold_error(scenario_pk$auc_to_last_obs, obs_auc_to_last_obs),
    3
  ),
  stringsAsFactors = FALSE
)

scenario_legend_labels <- paste0(
  error_summary$scenario,
  " (GMFE = ",
  round(error_summary$fold_error, 2),
  ")"
)

cat("\nScenario lookup:\n")
print(scenario_lookup, row.names = FALSE)

cat("\nScenario PK summary:\n")
print(
  transform(
    scenario_pk,
    auc_to_last_obs = round(auc_to_last_obs, 3)
  ),
  row.names = FALSE
)

cat("\nError summary:\n")
print(error_summary, row.names = FALSE)

# Plot specification ------------------------------------------------------

curve_spec <- data.frame(
  sim_path = rep(out_path, length(scenarios)),
  label = vapply(scenarios, \(x) x$id, character(1)),
  color = vapply(scenarios, \(x) x$color, character(1)),
  linewidth = vapply(scenarios, \(x) x$linewidth, numeric(1)),
  linetype = vapply(scenarios, \(x) x$linetype, character(1)),
  stringsAsFactors = FALSE
)

legend_spec <- build_legend_spec(
  curve_spec = curve_spec,
  include_observed = include_observed,
  obs_label = obs_label,
  obs_color = obs_color,
  obs_shape = obs_shape
)

p <- make_scenario_plot(
  plot_data = dc,
  legend_spec = legend_spec,
  title = plot_title,
  include_observed = include_observed,
  x_unit = ospUnits$Time$h,
  y_scale = y_scale,
  x_axis_label = x_axis_label,
  y_axis_label = y_axis_label,
  legend_position = legend_position
)

print(p)