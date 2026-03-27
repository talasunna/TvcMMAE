library(ospsuite)
library(ggplot2)

rm(list = ls())

#MMAE

# --- paths ---
sim_path <- "~/Desktop/Thesis/MMAE/Mouse1 MMAE i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|MMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.MMAE.IV__10mgKg_MMAE_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# --- helper functions ---
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

# --- scenario ---
scn <- list(
  name = "FU=0.83, Lipophilicity=2.855, SC=0.01",
  fu = 0.83,
  lipophilicity = 2.855,
  specific_clearance = 0.01
)

# =========================
# observed data extraction
# =========================
obs_df <- data.frame(
  time = obs_data$xValues,
  conc = obs_data$yValues
)

obs_df <- obs_df[is.finite(obs_df$time) & is.finite(obs_df$conc) & obs_df$conc > 0, ]

cat("Observed data set name:", obs_data$name, "\n")
cat("Observed x unit:", obs_data$xUnit, "\n")
cat("Observed y unit:", obs_data$yUnit, "\n\n")

# =========================
# simulation
# =========================
sim <- loadSimulation(sim_path)
setOutputs(out_path, sim)

apply_changes(
  sim,
  fu = scn$fu,
  lipophilicity = scn$lipophilicity,
  specific_clearance = scn$specific_clearance
)

res <- runSimulations(sim)[[1]]

# extract simulated time-concentration data
sim_out <- getOutputValues(
  simulationResults = res,
  quantitiesOrPaths = out_path,
  individualIds = 0
)

sim_df <- data.frame(
  time = sim_out$data$Time,
  conc = sim_out$data[[out_path]]
)

sim_df <- sim_df[is.finite(sim_df$time) & is.finite(sim_df$conc) & sim_df$conc > 0, ]

cat("Simulation output metadata:\n")
print(sim_out$metaData)

# =========================================
# IMPORTANT: make sure time units match
# getOutputValues() returns time in minutes
# if observed data are in hours, convert sim time
# =========================================
if (as.character(obs_data$xUnit) == "h") {
  sim_df$time <- sim_df$time / 60
}

# =========================================
# IMPORTANT: concentrations must also be in
# the same units before comparing Cmax/AUC
# inspect obs_data$yUnit and sim_out$metaData
# =========================================

# --- compute simple PK metrics ---
obs_cmax <- max(obs_df$conc)
sim_cmax <- max(sim_df$conc)

obs_auc_last <- auc_last(obs_df$time, obs_df$conc)
sim_auc_last <- auc_last(sim_df$time, sim_df$conc)

comparison <- data.frame(
  parameter = c("Cmax", "AUC_last"),
  observed = c(obs_cmax, obs_auc_last),
  predicted = c(sim_cmax, sim_auc_last)
)

comparison$fold_error <- ifelse(
  comparison$predicted >= comparison$observed,
  comparison$predicted / comparison$observed,
  comparison$observed / comparison$predicted
)

overall_gmfe <- gmfe(comparison$predicted, comparison$observed)

cat("\n====================\n")
cat("Scenario:", scn$name, "\n")
cat("====================\n\n")

print(comparison)
cat("\nOverall GMFE =", overall_gmfe, "\n")


dc <- DataCombined$new()
dc$addDataSets(obs_data)
dc$addSimulationResults(
  simulationResults = res,
  quantitiesOrPaths = out_path,
  names = scn$name
)

# --- plot in hours ---
plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)