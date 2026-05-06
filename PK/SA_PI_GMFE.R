library(ospsuite)
library(ospsuite.reportingengine)
library(ggplot2)

rm(list = ls())

options(width = 140)
options(ospsuite.plots.watermarkEnabled = FALSE)

# Paths -------------------------------------------------------------------

sim_path <- "~/Desktop/Thesis/PD/TGI.pkml"
out_path <- "Organism|Tumor|MMAE|totalMMAEtumor"
dataset_path <- "~/Desktop/Thesis/Data/Observed_data/tumorPK/Chang et al.totalMMAE.tumor.pkml"

# Plot options ------------------------------------------------------------

include_observed <- TRUE
plot_title <- "Plasma total mAb Concentration-Time Profiles"
x_axis_label <- "Time [h]"
y_axis_label <- "Concentration [\u03bcmol/l]"
legend_position <- c(0.8, 0.8)
legend.background <- element_rect(fill = "white", color = "grey70")
legend.key <- element_rect(fill = "white", color = NA)
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
    fu = NULL,
    lipophilicity = NULL,
    specific_clearance = NULL) {
  
  if (!is.null(tmmae_her2_kd)) {
    her2_kd <- getParameter("TvcMMAE-HER2-Chang et al|Kd", sim)
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
  
  if (!is.null(fu)) {
    fu_plasma <- getParameter("TvcMMAE|Fraction unbound (plasma, reference value)", sim)
    setParameterValues(fu_plasma, fu)
  }
  
  if (!is.null(lipophilicity)) {
    lip <- getParameter("TvcMMAE|Lipophilicity", sim)
    setParameterValues(lip, lipophilicity)
  }
  
  if (!is.null(specific_clearance)) {
    sc <- getParameter("TvcMMAE-Total Hepatic Clearance-Mouse|Specific clearance", sim)
    setParameterValues(sc, specific_clearance)
  }
  
  sim$solver$absTol <- 1e-12
  sim$solver$relTol <- 1e-8
  
  invisible(sim)
}

get_last_obs_time_h <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  max(obs_df$xValues, na.rm = TRUE)
}

calc_auc_to_time <- function(time_min, conc, end_time_min) {
  keep <- is.finite(time_min) &
    is.finite(conc) &
    time_min <= end_time_min
  
  time_min <- time_min[keep]
  conc <- conc[keep]
  
  ord <- order(time_min)
  time_min <- time_min[ord]
  conc <- conc[ord]
  
  if (length(time_min) < 2) {
    stop("Not enough points to calculate AUC.")
  }
  
  sum(diff(time_min) * (head(conc, -1) + tail(conc, -1)) / 2)
}

calculate_scenario_gmfe <- function(pred, obs) {
  pred <- as.numeric(pred)
  obs <- as.numeric(obs)
  
  if (length(pred) == 1 && length(obs) > 1) {
    pred <- rep(pred, length(obs))
  } else if (length(obs) == 1 && length(pred) > 1) {
    obs <- rep(obs, length(pred))
  } else if (length(pred) != length(obs)) {
    stop("`pred` and `obs` must have the same length, or one must have length 1.")
  }
  
  out <- rep(NA_real_, length(pred))
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  
  out[ok] <- vapply(
    seq_along(pred[ok]),
    \(i) {
      ospsuite.reportingengine::calculateGMFE(
        x = pred[ok][i],
        y = obs[ok][i]
      )
    },
    numeric(1)
  )
  
  out
}

build_legend_spec <- function(
    curve_spec,
    include_observed = FALSE,
    obs_label = "Observed data",
    obs_color = "orange",
    obs_shape = 16
) {
  color_values <- stats::setNames(curve_spec$color, curve_spec$name)
  linetype_values <- stats::setNames(curve_spec$linetype, curve_spec$name)
  linewidth_values <- stats::setNames(curve_spec$linewidth, curve_spec$name)
  shape_values <- stats::setNames(rep(NA, nrow(curve_spec)), curve_spec$name)
  
  label_lookup <- stats::setNames(curve_spec$label, curve_spec$name)
  
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
    label_lookup <- c(
      label_lookup,
      stats::setNames(obs_label, obs_label)
    )
    legend_order <- c(curve_spec$name, obs_label)
  } else {
    legend_order <- curve_spec$name
  }
  
  list(
    order = legend_order,
    labels = unname(label_lookup[legend_order]),
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
    legend_position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "grey70"),
    legend.key = element_rect(fill = "white", color = NA)
) {
  base_plot <- ospsuite::plotTimeProfile(
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
      breaks = legend_spec$order,
      labels = legend_spec$labels
    ) +
    scale_linetype_manual(
      values = legend_spec$linetype_values,
      breaks = legend_spec$order,
      labels = legend_spec$labels
    ) +
    scale_linewidth_manual(
      values = legend_spec$linewidth_values,
      breaks = legend_spec$order,
      labels = legend_spec$labels
    ) +
    scale_shape_manual(
      values = legend_spec$shape_values,
      breaks = legend_spec$order,
      labels = legend_spec$labels
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
      legend.title = element_blank(),
      legend.background = legend.background,
      legend.key = legend.key
    )
}

# Baseline sensitivity analysis -------------------------------------------

sim_sa <- loadSimulation(sim_path)
setOutputs(out_path, sim_sa)

sa_parameter_paths <- c(
  "**|TvcMMAE|Kd (FcRn) in Endosomal Space",
  "TvcMMAE-HER2-Chang et al|Kd",
  "HER2|Reference concentration",
  "TvcMMAE|Fraction unbound (plasma, reference value)",
  "TvcMMAE|Lipophilicity",
  "TvcMMAE-Total Hepatic Clearance-Mouse|Specific clearance"
  
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
    id = "Fraction unboud: 0.2",
    label = "FcRnKd: 1.07 \u03bcmol/l, HER2Kd: 1 \u03bcmol/l, HER2 ref: 0.7 \u03bcmol/l",
    fcrn_endo_kd = NULL,
    tmmae_her2_kd = NULL,
    herref = NULL,
    dose = NULL,
    fu = 0.2,
    lipophilicity = NULL,
    specific_clearance = NULL,
    color = "darkgreen",
    linewidth = 0.9,
    linetype = "solid"
  ),
  
  list(
    id = "Fraction unboud: 0.4",
    label = "FcRnKd: 1.07 \u03bcmol/l, HER2Kd: 8 \u03bcmol/l, HER2 ref: 0.25 \u03bcmol/l, fu: 0.5",
    fcrn_endo_kd = NULL,
    tmmae_her2_kd = NULL,
    herref = NULL,
    dose = NULL,
    fu = 0.4,
    lipophilicity = NULL,
    specific_clearance = NULL,
    color = "blue",
    linewidth = 0.9,
    linetype = "solid"
  ),
  
  list(
    id = "Fraction unboud: 0.8",
    label = "1",
    fcrn_endo_kd = NULL,
    tmmae_her2_kd = NULL,
    herref = NULL,
    dose = NULL,
    fu = 0.8,
    lipophilicity = NULL,
    specific_clearance = NULL,
    color = "red",
    linewidth = 0.9,
    linetype = "solid"
  )
)

# Observed reference metric -----------------------------------------------

last_obs_time_h <- get_last_obs_time_h(obs_data)
auc_end_time_min <- last_obs_time_h * 60

obs_df <- dataSetToDataFrame(obs_data)

obs_auc_to_last_obs <- calc_auc_to_time(
  time_min = obs_df$xValues * 60,
  conc = obs_df$yValues,
  end_time_min = auc_end_time_min
)

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
    fu = scn$fu,
    lipophilicity = scn$lipophilicity,
    specific_clearance = scn$specific_clearance
  )
  
  res <- runSimulations(sim)[[1]]
  
  sim_values <- getOutputValues(
    simulationResults = res,
    quantitiesOrPaths = out_path
  )
  
  sim_df <- sim_values$data
  
  sim_auc_to_last_obs <- calc_auc_to_time(
    time_min = sim_df$Time,
    conc = sim_df[[out_path]],
    end_time_min = auc_end_time_min
  )
  
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

scenario_gmfe <- calculate_scenario_gmfe(
  pred = scenario_pk$auc_to_last_obs,
  obs = obs_auc_to_last_obs
)

error_summary <- data.frame(
  scenario = scenario_pk$scenario,
  observed_auc_to_last_obs = round(
    rep(obs_auc_to_last_obs, nrow(scenario_pk)),
    3
  ),
  predicted_auc_to_last_obs = round(scenario_pk$auc_to_last_obs, 3),
  gmfe = round(scenario_gmfe, 3),
  stringsAsFactors = FALSE
)

scenario_legend_labels <- paste0(
  error_summary$scenario,
  " (GMFE: ",
  formatC(error_summary$gmfe, format = "f", digits = 2),
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
  name = vapply(scenarios, \(x) x$id, character(1)),
  label = scenario_legend_labels,
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
  legend_position = legend_position,
  legend.background = legend.background,
  legend.key = legend.key
)

print(p)