library(ospsuite)
library(ggplot2)

rm(list = ls())

##  TMMAE
# --- paths ---
sim_path <- "~/Desktop/Thesis/TMMAE/Mouse1 _DC_i.v. _10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/TMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# optional: inspect what was loaded
print(obs_data)

# --------------------------------------------------
# optional MW only needed if concentration units are
# mass-based (e.g. mg/L, ug/mL) instead of molar
# --------------------------------------------------
manual_mw_g_per_mol <- 155000

# -----------------------------
# helper functions
# -----------------------------
gmfe <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & pred > 0 & obs > 0
  if (!any(ok)) stop("No valid positive observed/predicted pairs for GMFE.")
  10^(mean(abs(log10(pred[ok] / obs[ok]))))
}

auc_last <- function(time, conc) {
  o <- order(time)
  time <- time[o]
  conc <- conc[o]
  sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)
}

get_pk_value <- function(pk_df, quantity_path, parameter_name) {
  row <- pk_df[pk_df$QuantityPath == quantity_path & pk_df$Parameter == parameter_name, ]
  if (nrow(row) == 0) stop(paste("PK parameter not found:", parameter_name))
  row$Value[1]
}

get_pk_unit <- function(pk_df, quantity_path, parameter_name) {
  row <- pk_df[pk_df$QuantityPath == quantity_path & pk_df$Parameter == parameter_name, ]
  if (nrow(row) == 0) stop(paste("PK parameter not found:", parameter_name))
  as.character(row$Unit[1])
}

normalize_unit <- function(u) {
  u <- as.character(u)
  u <- gsub("µ", "u", u, fixed = TRUE)
  u <- tolower(trimws(u))
  u
}

get_mw <- function(obs_data, manual_mw_g_per_mol = 155000) {
  mw <- suppressWarnings(as.numeric(obs_data$molecularWeight))
  if (!is.finite(mw)) {
    mw <- manual_mw_g_per_mol
  }
  mw
}

convert_time_to_h <- function(value, unit) {
  u <- normalize_unit(unit)
  
  if (u %in% c("h", "hr", "hrs", "hour", "hours")) return(value)
  if (u %in% c("min", "mins", "minute", "minutes")) return(value / 60)
  if (u %in% c("s", "sec", "secs", "second", "seconds")) return(value / 3600)
  
  warning(paste("Unexpected time unit:", unit, "- leaving unchanged"))
  value
}

convert_conc_to_umol_per_l <- function(value, unit, mw_g_per_mol = NA_real_) {
  u <- normalize_unit(unit)
  
  if (u == "umol/l") return(value)
  if (u == "nmol/l") return(value / 1000)
  if (u == "mmol/l") return(value * 1000)
  
  # mass-based units need MW
  if (!is.finite(mw_g_per_mol)) {
    stop(paste(
      "Concentration unit", unit,
      "requires molecular weight conversion. Set manual_mw_g_per_mol."
    ))
  }
  
  if (u %in% c("mg/l", "ug/ml")) return(value * 1000 / mw_g_per_mol)
  if (u %in% c("ug/l", "ng/ml")) return(value / mw_g_per_mol)
  
  warning(paste("Unexpected concentration unit:", unit, "- leaving unchanged"))
  value
}

convert_auc_to_umol_h_per_l <- function(value, unit, mw_g_per_mol = NA_real_) {
  u <- normalize_unit(unit)
  
  if (u == "umol*h/l") return(value)
  if (u == "umol*min/l") return(value / 60)
  
  if (u == "nmol*h/l") return(value / 1000)
  if (u == "nmol*min/l") return(value / 1000 / 60)
  
  if (u == "mmol*h/l") return(value * 1000)
  if (u == "mmol*min/l") return(value * 1000 / 60)
  
  # mass-based units need MW
  if (!is.finite(mw_g_per_mol)) {
    stop(paste(
      "AUC unit", unit,
      "requires molecular weight conversion. Set manual_mw_g_per_mol."
    ))
  }
  
  if (u %in% c("mg*h/l", "ug*h/ml")) return(value * 1000 / mw_g_per_mol)
  if (u %in% c("mg*min/l", "ug*min/ml")) return(value * 1000 / mw_g_per_mol / 60)
  
  if (u %in% c("ug*h/l", "ng*h/ml")) return(value / mw_g_per_mol)
  if (u %in% c("ug*min/l", "ng*min/ml")) return(value / mw_g_per_mol / 60)
  
  warning(paste("Unexpected AUC unit:", unit, "- leaving unchanged"))
  value
}

# --- function that applies one scenario to a simulation ---
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

# --- scenarios ---
scenarios <- list(
  list(
    name = "FcRnKd: 1.07 μmol/l, HER2Kd: 3.92 μmol/l, HER2 ref: 1 μmol/l, Dose: 15 mg/kg, FU = 0.5",
    fcrn_endo_kd = 1,
    tmmae_her2_kd = 3.92,
    herref = 1,
    dose = NULL,
    fu = 0.2
  )
)

# -----------------------------
# observed data extraction
# -----------------------------
mw_g_per_mol <- get_mw(obs_data, manual_mw_g_per_mol)

obs_time_h <- convert_time_to_h(obs_data$xValues, obs_data$xUnit)

obs_df <- data.frame(
  time_h = obs_time_h,
  conc_raw = obs_data$yValues
)

obs_df <- obs_df[
  is.finite(obs_df$time_h) &
    is.finite(obs_df$conc_raw) &
    obs_df$conc_raw > 0,
]

obs_cmax_raw <- max(obs_df$conc_raw)
obs_tmax_h <- obs_df$time_h[which.max(obs_df$conc_raw)[1]]
obs_auc_last_raw_h <- auc_last(obs_df$time_h, obs_df$conc_raw)

obs_cmax_umol_per_l <- convert_conc_to_umol_per_l(
  obs_cmax_raw,
  obs_data$yUnit,
  mw_g_per_mol
)

# same concentration conversion factor applies to AUC because time is already in h
obs_auc_last_umol_h_per_l <- convert_conc_to_umol_per_l(
  obs_auc_last_raw_h,
  obs_data$yUnit,
  mw_g_per_mol
)

cat("Observed data set name:", obs_data$name, "\n")
cat("Observed x unit:", obs_data$xUnit, "\n")
cat("Observed y unit:", obs_data$yUnit, "\n")
cat("Observed MW from dataset:", obs_data$molecularWeight, "\n")
cat("MW used for conversion:", mw_g_per_mol, "\n\n")

cat("Observed PK:\n")
cat("Cmax =", obs_cmax_umol_per_l, "umol/L\n")
cat("Tmax =", obs_tmax_h, "h\n")
cat("AUC_last =", obs_auc_last_umol_h_per_l, "umol*h/L\n\n")

# -----------------------------
# add user-defined AUC over observed time window
# -----------------------------
removeAllUserDefinedPKParameters()

my_auc_last <- addUserDefinedPKParameter(
  name = "AUC_last_obsWindow",
  standardPKParameter = StandardPKParameter$AUC_tEnd
)

# OSP expects minutes here
my_auc_last$startTime <- 0
my_auc_last$endTime <- max(obs_df$time_h) * 60

# --- combine observed + simulated data ---
dc <- DataCombined$new()
dc$addDataSets(obs_data)

gmfe_results <- list()

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
  
  dc$addSimulationResults(
    simulationResults = res,
    quantitiesOrPaths = out_path,
    names = scn$name
  )
  
  # -----------------------------
  # PK analysis for this scenario
  # -----------------------------
  pk_analysis <- calculatePKAnalyses(results = res)
  pk_df <- pkAnalysesToDataFrame(pk_analysis)
  
  sim_cmax <- get_pk_value(pk_df, out_path, "C_max")
  sim_cmax_unit <- get_pk_unit(pk_df, out_path, "C_max")
  
  sim_tmax <- get_pk_value(pk_df, out_path, "t_max")
  sim_tmax_unit <- get_pk_unit(pk_df, out_path, "t_max")
  
  sim_auc_last <- get_pk_value(pk_df, out_path, "AUC_last_obsWindow")
  sim_auc_last_unit <- get_pk_unit(pk_df, out_path, "AUC_last_obsWindow")
  
  sim_cmax_umol_per_l <- convert_conc_to_umol_per_l(
    sim_cmax,
    sim_cmax_unit,
    mw_g_per_mol
  )
  
  sim_tmax_h <- convert_time_to_h(
    sim_tmax,
    sim_tmax_unit
  )
  
  sim_auc_last_umol_h_per_l <- convert_auc_to_umol_h_per_l(
    sim_auc_last,
    sim_auc_last_unit,
    mw_g_per_mol
  )
  
  comparison <- data.frame(
    parameter = c("Cmax", "AUC_last"),
    observed = c(
      obs_cmax_umol_per_l,
      obs_auc_last_umol_h_per_l
    ),
    predicted = c(
      sim_cmax_umol_per_l,
      sim_auc_last_umol_h_per_l
    )
  )
  
  comparison$fold_error <- ifelse(
    comparison$predicted >= comparison$observed,
    comparison$predicted / comparison$observed,
    comparison$observed / comparison$predicted
  )
  
  overall_gmfe <- gmfe(comparison$predicted, comparison$observed)
  
  tmax_table <- data.frame(
    parameter = "Tmax",
    observed_h = obs_tmax_h,
    predicted_h = sim_tmax_h,
    fold_error = ifelse(
      sim_tmax_h >= obs_tmax_h,
      sim_tmax_h / obs_tmax_h,
      obs_tmax_h / sim_tmax_h
    )
  )
  
  cat("\n====================\n")
  cat("Scenario:", scn$name, "\n")
  cat("====================\n\n")
  
  cat("Simulated PK from OSP Suite:\n")
  cat("C_max =", sim_cmax_umol_per_l, "umol/L", "(raw:", sim_cmax, sim_cmax_unit, ")\n")
  cat("t_max =", sim_tmax_h, "h", "(raw:", sim_tmax, sim_tmax_unit, ")\n")
  cat("AUC_last_obsWindow =", sim_auc_last_umol_h_per_l, "umol*h/L",
      "(raw:", sim_auc_last, sim_auc_last_unit, ")\n\n")
  
  cat("Comparison table:\n")
  print(comparison)
  
  cat("\nTmax comparison:\n")
  print(tmax_table)
  
  cat("\nOverall PK GMFE (Cmax + AUC_last) =", overall_gmfe, "\n")
  
  gmfe_results[[i]] <- data.frame(
    scenario = scn$name,
    obs_cmax_umol_per_l = obs_cmax_umol_per_l,
    sim_cmax_umol_per_l = sim_cmax_umol_per_l,
    obs_tmax_h = obs_tmax_h,
    sim_tmax_h = sim_tmax_h,
    obs_auc_last_umol_h_per_l = obs_auc_last_umol_h_per_l,
    sim_auc_last_umol_h_per_l = sim_auc_last_umol_h_per_l,
    gmfe_pk = overall_gmfe,
    tmax_fold_error = tmax_table$fold_error
  )
}

gmfe_results_df <- do.call(rbind, gmfe_results)

cat("\n========================================\n")
cat("Summary across scenarios\n")
cat("========================================\n")
print(gmfe_results_df)

# --- plot in hours ---
plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)