library(ospsuite)
library(ggplot2)

rm(list = ls())

# --- paths ---
sim_path <- "~/Desktop/Thesis/TMMAE/Mouse1 ADC i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/TMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)
print(obs_data)

# =========================
# baseline sensitivity analysis
# =========================

sim_sa <- loadSimulation(sim_path)
setOutputs(out_path, sim_sa)
#potentialSAParameters <- potentialVariableParameterPathsFor(sim_sa)
#print(potentialSAParameters)

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

sens_cmax <- sa_result$allPKParameterSensitivitiesFor(
  pkParameterName = "C_max",
  outputPath = out_path,
  totalSensitivityThreshold = 1
)

sens_auc <- sa_result$allPKParameterSensitivitiesFor(
  pkParameterName = "AUC_tEnd",
  outputPath = out_path,
  totalSensitivityThreshold = 1
)

sens_cl <- sa_result$allPKParameterSensitivitiesFor(
  pkParameterName = "CL",
  outputPath = out_path,
  totalSensitivityThreshold = 1
)

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

sa_df <- rbind(
  to_df(sens_cmax, "C_max"),
  to_df(sens_auc, "AUC_tEnd"),
  to_df(sens_cl, "CL")
)

sa_df <- sa_df[order(sa_df$PK_Parameter, -abs(sa_df$Sensitivity)), ]
print(sa_df)

# =========================
# helper functions
# =========================

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

gmfe <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  if (!any(ok)) stop("No valid positive observed/predicted pairs for GMFE.")
  10^(mean(abs(log10(pred[ok] / obs[ok]))))
}

get_last_obs_time <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  max(obs_df$xValues, na.rm = TRUE)
}

get_obs_cmax <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  max(obs_df$yValues, na.rm = TRUE)
}

get_obs_auc_last_obs_window <- function(obs_data) {
  obs_df <- dataSetToDataFrame(obs_data)
  obs_df <- obs_df[order(obs_df$xValues), ]
  
  x <- obs_df$xValues
  y <- obs_df$yValues
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  # trapezoidal AUC in µmol*h/l
  auc_h <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  
  # convert to µmol*min/l to match OSP PK analysis output
  auc_min <- auc_h * 60
  
  auc_min
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

get_sim_auc_last_obs_window <- function(sim_results, output_path, end_time_h) {
  removeAllUserDefinedPKParameters()
  
  auc_obs_window <- addUserDefinedPKParameter(
    name = "AUC_last_obsWindow",
    standardPKParameter = StandardPKParameter$AUC_tEnd
  )
  
  auc_obs_window$startTime <- 0
  auc_obs_window$endTime <- end_time_h *60
  
  pk <- calculatePKAnalyses(results = sim_results)
  pk_df <- pkAnalysesToDataFrame(pk)
  
  row <- pk_df[pk_df$QuantityPath == output_path & pk_df$Parameter == "AUC_last_obsWindow", ]
  
  if (nrow(row) == 0) {
    stop("User-defined parameter AUC_last_obsWindow was not found in PK analysis output.")
  }
  
  row$Value[[1]]
}
# =========================
# scenarios
# =========================

scenarios <- list(
  list(
    name = "FcRnKd: 1.07 μmol/l, HER2Kd: 3.92 μmol/l, HER2 ref: 1 μmol/l",
    fcrn_endo_kd = 1,
    tmmae_her2_kd = 3.92,
    herref = 1,
    dose = NULL,
    fu = NULL
  )
  # list(name = "FU = 0.5", fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL, fu = 0.5),
  # list(name = "FU = 0.8", fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL, fu = 0.8)
)

# ==========================
# observed reference metrics
# ==========================

last_obs_time_h <- get_last_obs_time(obs_data)
obs_cmax <- get_obs_cmax(obs_data)
obs_auc_last_obs_window <- get_obs_auc_last_obs_window(obs_data)

cat("Last observed time (h):", last_obs_time_h, "\n")
cat("Observed Cmax:", obs_cmax, "\n")
cat("Observed AUC_last_obsWindow (µmol*min/l):", obs_auc_last_obs_window, "\n")

# =================================
# combined plot + scenario PK table
# =================================

dc <- DataCombined$new()
dc$addDataSets(obs_data)

scenario_pk <- data.frame(
  Scenario = character(),
  C_max = numeric(),
  AUC_last_obsWindow = numeric(),
  stringsAsFactors = FALSE
)

for (scn in scenarios) {
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
  sim_vals <- getOutputValues(res, quantitiesOrPaths = out_path)
  cat("Max simulated time stored (min):", max(sim_vals$data$Time, na.rm = TRUE), "\n")
  cat("Requested AUC end time (min):", last_obs_time_h * 60, "\n")
  
  sim_cmax <- get_sim_pk_value(res, out_path, "C_max")
  sim_auc_last_obs_window <- get_sim_auc_last_obs_window(res, out_path, last_obs_time_h)
  
  scenario_pk <- rbind(
    scenario_pk,
    data.frame(
      Scenario = scn$name,
      C_max = sim_cmax,
      AUC_last_obsWindow = sim_auc_last_obs_window,
      stringsAsFactors = FALSE
    )
  )
  
  dc$addSimulationResults(
    simulationResults = res,
    quantitiesOrPaths = out_path,
    names = scn$name
  )
}

print(scenario_pk)

# =========================
# GMFE
# =========================

gmfe <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  if (!any(ok)) stop("No valid positive observed/predicted pairs for GMFE.")
  10^(mean(abs(log10(pred[ok] / obs[ok]))))
}

fold_error <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  out <- rep(NA_real_, length(pred))
  out[ok] <- 10^(abs(log10(pred[ok] / obs[ok])))
  out
}

gmfe_tables <- list()

for (i in seq_len(nrow(scenario_pk))) {
  
  pred_cmax <- scenario_pk$C_max[i]
  pred_auc  <- scenario_pk$AUC_last_obsWindow[i]
  
  result_table <- data.frame(
    Parameter = c("C_max", "AUC_last_obsWindow"),
    Observed = c(obs_cmax, obs_auc_last_obs_window),
    Predicted = c(pred_cmax, pred_auc),
    Fold_Error = c(
      fold_error(pred_cmax, obs_cmax),
      fold_error(pred_auc, obs_auc_last_obs_window)
    ),
    stringsAsFactors = FALSE
  )
  
  overall_gmfe <- gmfe(
    pred = result_table$Predicted,
    obs  = result_table$Observed
  )
  
  gmfe_tables[[scenario_pk$Scenario[i]]] <- list(
    table = result_table,
    overall_gmfe = overall_gmfe
  )
}

for (scn_name in names(gmfe_tables)) {
  cat("\n=============================\n")
  cat("Scenario:", scn_name, "\n")
  cat("=============================\n")
  print(gmfe_tables[[scn_name]]$table, row.names = FALSE)
  cat("\nOverall GMFE:", gmfe_tables[[scn_name]]$overall_gmfe, "\n")
}

# =========================
# plot
# =========================

plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)