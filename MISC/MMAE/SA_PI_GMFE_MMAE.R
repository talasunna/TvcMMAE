library(ospsuite)

rm(list = ls())
options(width = 140)

# --- paths ---
sim_path <- "~/Desktop/Thesis/MMAE/Mouse1 MMAE i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|MMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.MMAE.IV__10mgKg_MMAE_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)
print(obs_data)

# =========================
# helper functions
# =========================

to_df <- function(x, pk_name) {
  if (length(x) == 0) {
    return(data.frame(
      PK_Parameter = character(0),
      Parameter = character(0),
      Sensitivity = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    PK_Parameter = rep(pk_name, length(x)),
    Parameter = sapply(x, function(y) y$parameterName),
    Sensitivity = sapply(x, function(y) y$value),
    stringsAsFactors = FALSE
  )
}

apply_changes <- function(sim,
                          fu = NULL,
                          lipophilicity = NULL,
                          specific_clearance = NULL) {
  
  if (!is.null(fu)) {
    fu_plasma <- getParameter("MMAE|Fraction unbound (plasma, reference value)", sim)
    setParameterValues(fu_plasma, fu)
  }
  
  if (!is.null(lipophilicity)) {
    lip <- getParameter("MMAE|Lipophilicity", sim)
    setParameterValues(lip, lipophilicity)
  }
  
  if (!is.null(specific_clearance)) {
    sc <- getParameter("MMAE-Total Hepatic Clearance-Mouse|Specific clearance", sim)
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

get_obs_auc_to_last_obs <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  obs_df <- obs_df[order(obs_df$xValues), ]
  
  x <- obs_df$xValues
  y <- obs_df$yValues
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  # Assumes xValues are in hours and yValues are in µmol/l
  # Trapezoidal AUC first in µmol*h/l, then converted to µmol*min/l
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
  
  row <- pk_df[pk_df$QuantityPath == output_path & pk_df$Parameter == parameter_name, ]
  
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

# =========================
# baseline sensitivity analysis
# =========================
# Sensitivity analysis uses standard PK parameters available in the SA output.

sim_sa <- loadSimulation(sim_path)
setOutputs(out_path, sim_sa)

sa_parameter_paths <- c(
  "MMAE|Fraction unbound (plasma, reference value)",
  "MMAE|Lipophilicity",
  "MMAE-Total Hepatic Clearance-Mouse|Specific clearance"
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

sa_df <- sa_df[order(sa_df$PK_Parameter, -abs(sa_df$Sensitivity)), ]

cat("\nSensitivity summary:\n")
print(sa_df, row.names = FALSE)

# =========================
# scenarios
# =========================

scenarios <- list(
  list(
    id = "S1",
    label = "FU = 0.83, Lip = 3, SC = 0.006 ",
    fu = 0.83,
    lipophilicity = 3,
    specific_clearance = 0.0255
  )
  
  # ,
  # list(
  #   id = "S2",
  #   label = "FU = 0.5, Lip = 3.0, SC = 0.0255",
  #   fu = 0.5,
  #   lipophilicity = 3.0,
  #   specific_clearance = 0.015
  # ),
  #
  # list(
  #   id = "S3",
  #   label = "FU = 0.7, Lip = 2.5, SC = 0.03",
  #   fu = 0.7,
  #   lipophilicity = 2.5,
  #   specific_clearance = 0.03
  # )
)

# ==========================
# observed reference metric
# ==========================

last_obs_time_h <- get_last_obs_time_h(obs_data)
obs_auc_to_last_obs <- get_obs_auc_to_last_obs(obs_data)

define_auc_to_last_obs_parameter(last_obs_time_h)

cat("\nObserved reference:\n")
cat("Last observed time (h):", round(last_obs_time_h, 3), "\n")
cat("Observed AUC_to_lastObs (µmol*min/l):", round(obs_auc_to_last_obs, 3), "\n")

# =================================
# combined plot + scenario PK table
# =================================

dc <- DataCombined$new()
dc$addDataSets(obs_data)

scenario_pk_list <- vector("list", length(scenarios))
scenario_lookup_list <- vector("list", length(scenarios))

for (i in seq_along(scenarios)) {
  scn <- scenarios[[i]]
  
  sim <- loadSimulation(sim_path)
  setOutputs(out_path, sim)
  
  apply_changes(
    sim,
    fu = scn$fu,
    lipophilicity = scn$lipophilicity,
    specific_clearance = scn$specific_clearance
  )
  
  res <- runSimulations(sim)[[1]]
  
  sim_auc_to_last_obs <- get_sim_pk_value(res, out_path, "AUC_to_lastObs")
  
  scenario_pk_list[[i]] <- data.frame(
    Scenario = scn$id,
    AUC_to_lastObs = sim_auc_to_last_obs,
    stringsAsFactors = FALSE
  )
  
  scenario_lookup_list[[i]] <- data.frame(
    Scenario = scn$id,
    Description = scn$label,
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

# =========================
# error analysis
# =========================

error_summary <- data.frame(
  Scenario = scenario_pk$Scenario,
  Observed_AUC_to_lastObs = round(rep(obs_auc_to_last_obs, nrow(scenario_pk)), 3),
  Predicted_AUC_to_lastObs = round(scenario_pk$AUC_to_lastObs, 3),
  Fold_Error = round(fold_error(scenario_pk$AUC_to_lastObs, obs_auc_to_last_obs), 3),
  stringsAsFactors = FALSE
)

cat("\nScenario lookup:\n")
print(scenario_lookup, row.names = FALSE)

cat("\nScenario PK summary:\n")
print(
  transform(
    scenario_pk,
    AUC_to_lastObs = round(AUC_to_lastObs, 3)
  ),
  row.names = FALSE
)

cat("\nError summary:\n")
print(error_summary, row.names = FALSE)

# =========================
# plot
# =========================

plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)