library(ospsuite)

rm(list = ls())
options(width = 140)

# --- paths ---
sim_path <- "~/Desktop/Thesis/TvcMMAE/Mouse1 ADC i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TvcMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/TvcMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

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
                          fcrn_endo_kd = NULL,
                          tmmae_her2_kd = NULL,
                          herref = NULL,
                          dose = NULL,
                          fu = NULL) {
  
  if (!is.null(fu)) {
    fu_plasma <- getParameter("TMMAE|Fraction unbound (plasma, reference value)", sim)
    setParameterValues(fu_plasma, fu)
  }
  
  if (!is.null(tmmae_her2_kd)) {
    her2_kd <- getParameter("TMMAE-HER2-Chang et al|Kd", sim)
    setParameterValues(her2_kd, tmmae_her2_kd)
  }
  
  if (!is.null(fcrn_endo_kd)) {
    fcrnkd <- getParameter("**|TMMAE|Kd (FcRn) in Endosomal Space", sim)
    setParameterValues(fcrnkd, fcrn_endo_kd)
  }
  
  if (!is.null(dose)) {
    dose_par <- getParameter("Events|10 mg/kg|Application_1|ProtocolSchemaItem|DosePerBodyWeight", sim)
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
  "**|TMMAE|Kd (FcRn) in Endosomal Space",
  "TMMAE-HER2-Chang et al|Kd",
  "HER2|Reference concentration",
  "TMMAE|Fraction unbound (plasma, reference value)"
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
    label = "FcRnKd: 1.07 µmol/l, HER2Kd: 1 µmol/l, HER2 ref: 0.7 µmol/l",
    fcrn_endo_kd = 1.07,
    tmmae_her2_kd = 1,
    herref = 0.7,
    dose = NULL,
    fu = NULL
  )
  
  # ,
  # list(
  #   id = "S2",
  #   label = "FcRnKd: 1.07 µmol/l, HER2Kd: 8 µmol/l, HER2 ref: 0.25 µmol/l, fu: 0.5",
  #   fcrn_endo_kd = 1.07,
  #   tmmae_her2_kd = 8,
  #   herref = 0.25,
  #   dose = NULL,
  #   fu = 0.5
  # ),
  #
  # list(
  #   id = "S3",
  #   label = "FcRnKd: 1.07 µmol/l, HER2Kd: 10 µmol/l, HER2 ref: 0.1 µmol/l, fu: 0.75",
  #   fcrn_endo_kd = 1.07,
  #   tmmae_her2_kd = 10,
  #   herref = 0.1,
  #   dose = NULL,
  #   fu = 0.75
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
    fcrn_endo_kd  = scn$fcrn_endo_kd,
    tmmae_her2_kd = scn$tmmae_her2_kd,
    herref        = scn$herref,
    dose          = scn$dose,
    fu            = scn$fu
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