library(ospsuite)
library(ggplot2)

rm(list = ls())

# -----------------------------
# paths
# -----------------------------
sim_path <- "~/Desktop/Thesis/MMAE/Mouse1 MMAE i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|MMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.MMAE.IV__10mgKg_MMAE_Mouse.pkml"

# -----------------------------
# load observed data
# -----------------------------
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# -----------------------------
# USER INPUTS
# -----------------------------
dose_mg_per_kg <- 10
manual_mw_g_per_mol <- 717.98

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

estimate_lambda_z <- function(time, conc, min_points = 3) {
  o <- order(time)
  time <- time[o]
  conc <- conc[o]
  
  keep <- is.finite(time) & is.finite(conc) & conc > 0
  time <- time[keep]
  conc <- conc[keep]
  
  if (length(time) < min_points) {
    return(list(lambda_z = NA_real_, adj_r2 = NA_real_, used_idx = integer(0)))
  }
  
  idx_after_cmax <- which(seq_along(time) > which.max(conc))
  
  if (length(idx_after_cmax) < min_points) {
    return(list(lambda_z = NA_real_, adj_r2 = NA_real_, used_idx = integer(0)))
  }
  
  best_adj_r2 <- -Inf
  best_lambda <- NA_real_
  best_idx <- integer(0)
  
  for (k in min_points:length(idx_after_cmax)) {
    idx <- tail(idx_after_cmax, k)
    fit <- lm(log(conc[idx]) ~ time[idx])
    slope <- coef(fit)[2]
    adj_r2 <- summary(fit)$adj.r.squared
    
    if (is.finite(adj_r2) && is.finite(slope) && slope < 0 && adj_r2 > best_adj_r2) {
      best_adj_r2 <- adj_r2
      best_lambda <- -as.numeric(slope)
      best_idx <- idx
    }
  }
  
  list(lambda_z = best_lambda, adj_r2 = best_adj_r2, used_idx = best_idx)
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

convert_time_to_h <- function(value, unit) {
  u <- normalize_unit(unit)
  if (u == "h") return(value)
  if (u == "min") return(value / 60)
  warning(paste("Unexpected time unit:", unit, "- leaving unchanged"))
  value
}

convert_auc_to_umol_h_per_l <- function(value, unit) {
  u <- normalize_unit(unit)
  if (u %in% c("umol*h/l", "umol*h/l ")) return(value)
  if (u %in% c("umol*min/l", "umol*min/l ")) return(value / 60)
  warning(paste("Unexpected AUC unit:", unit, "- leaving unchanged"))
  value
}

convert_cl_to_l_per_min_per_kg <- function(value, unit) {
  u <- normalize_unit(unit)
  if (u == "l/min/kg") return(value)
  if (u == "ml/min/kg") return(value / 1000)
  warning(paste("Unexpected CL unit:", unit, "- leaving unchanged"))
  value
}

convert_vd_to_l_per_kg <- function(value, unit) {
  u <- normalize_unit(unit)
  if (u == "l/kg") return(value)
  if (u == "ml/kg") return(value / 1000)
  warning(paste("Unexpected Vd unit:", unit, "- leaving unchanged"))
  value
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

# -----------------------------
# scenario
# -----------------------------
scn <- list(
  name = "FU=0.83, Lip=2, SC=0.025",
  fu = 0.83,
  lipophilicity = 2,
  specific_clearance = 0.025
)

# -----------------------------
# observed data extraction
# -----------------------------
obs_df <- data.frame(
  time = obs_data$xValues,
  conc = obs_data$yValues
)

obs_df <- obs_df[
  is.finite(obs_df$time) &
    is.finite(obs_df$conc) &
    obs_df$conc > 0,
]

cat("Observed data set name:", obs_data$name, "\n")
cat("Observed x unit:", obs_data$xUnit, "\n")
cat("Observed y unit:", obs_data$yUnit, "\n")
cat("Observed MW from dataset:", obs_data$molecularWeight, "\n\n")

# molecular weight
mw_g_per_mol <- obs_data$molecularWeight
if (is.null(mw_g_per_mol) || !is.finite(mw_g_per_mol)) {
  mw_g_per_mol <- manual_mw_g_per_mol
}
if (!is.finite(mw_g_per_mol)) {
  stop("No molecular weight available. Set manual_mw_g_per_mol.")
}

# observed PK
obs_cmax <- max(obs_df$conc)
obs_tmax_h <- obs_df$time[which.max(obs_df$conc)[1]]
obs_auc_last_h <- auc_last(obs_df$time, obs_df$conc)

lz_obs <- estimate_lambda_z(obs_df$time, obs_df$conc)
obs_lambda_z_per_h <- lz_obs$lambda_z

if (!is.finite(obs_lambda_z_per_h) || obs_lambda_z_per_h <= 0) {
  stop("Could not estimate observed lambda_z reliably.")
}

obs_clast <- tail(obs_df$conc[order(obs_df$time)], 1)
obs_auc_inf_h <- obs_auc_last_h + obs_clast / obs_lambda_z_per_h

# dose mg/kg -> umol/kg
dose_umol_per_kg <- dose_mg_per_kg * 1000 / mw_g_per_mol

# observed CL and Vd
obs_cl_L_per_h_per_kg <- dose_umol_per_kg / obs_auc_inf_h
obs_cl_L_per_min_per_kg <- obs_cl_L_per_h_per_kg / 60
obs_vd_L_per_kg <- obs_cl_L_per_h_per_kg / obs_lambda_z_per_h

# -----------------------------
# simulation
# -----------------------------
sim <- loadSimulation(sim_path)
setOutputs(out_path, sim)

apply_changes(
  sim,
  fu = scn$fu,
  lipophilicity = scn$lipophilicity,
  specific_clearance = scn$specific_clearance
)

# user-defined PK parameter must exist before PK analysis
removeAllUserDefinedPKParameters()

my_auc_last <- addUserDefinedPKParameter(
  name = "AUC_last_obsWindow",
  standardPKParameter = StandardPKParameter$AUC_tEnd
)

# OSP expects minutes here
my_auc_last$startTime <- 0
my_auc_last$endTime <- max(obs_df$time) * 60

res <- runSimulations(sim)[[1]]

# -----------------------------
# OSP PK analysis
# -----------------------------
pk_analysis <- calculatePKAnalyses(results = res)
pk_df <- pkAnalysesToDataFrame(pk_analysis)

cat("Simulated PK parameters from OSP Suite:\n")
print(pk_df[pk_df$QuantityPath == out_path, ])

# raw simulation PK values
sim_cmax <- get_pk_value(pk_df, out_path, "C_max")
sim_cmax_unit <- get_pk_unit(pk_df, out_path, "C_max")

sim_tmax <- get_pk_value(pk_df, out_path, "t_max")
sim_tmax_unit <- get_pk_unit(pk_df, out_path, "t_max")

sim_auc_last <- get_pk_value(pk_df, out_path, "AUC_last_obsWindow")
sim_auc_last_unit <- get_pk_unit(pk_df, out_path, "AUC_last_obsWindow")

sim_auc_inf <- get_pk_value(pk_df, out_path, "AUC_inf")
sim_auc_inf_unit <- get_pk_unit(pk_df, out_path, "AUC_inf")

sim_cl <- get_pk_value(pk_df, out_path, "CL")
sim_cl_unit <- get_pk_unit(pk_df, out_path, "CL")

sim_vd <- get_pk_value(pk_df, out_path, "Vd")
sim_vd_unit <- get_pk_unit(pk_df, out_path, "Vd")

# -----------------------------
# convert simulation units
# -----------------------------
sim_tmax_h <- convert_time_to_h(sim_tmax, sim_tmax_unit)
sim_auc_last_h <- convert_auc_to_umol_h_per_l(sim_auc_last, sim_auc_last_unit)
sim_auc_inf_h <- convert_auc_to_umol_h_per_l(sim_auc_inf, sim_auc_inf_unit)
sim_cl_L_per_min_per_kg <- convert_cl_to_l_per_min_per_kg(sim_cl, sim_cl_unit)
sim_vd_L_per_kg <- convert_vd_to_l_per_kg(sim_vd, sim_vd_unit)

# -----------------------------
# comparison / GMFE
# -----------------------------
comparison <- data.frame(
  parameter = c("Cmax", "AUC_last", "CL", "Vd"),
  observed = c(
    obs_cmax,
    obs_auc_last_h,
    obs_cl_L_per_min_per_kg,
    obs_vd_L_per_kg
  ),
  predicted = c(
    sim_cmax,
    sim_auc_last_h,
    sim_cl_L_per_min_per_kg,
    sim_vd_L_per_kg
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

# -----------------------------
# print results
# -----------------------------
cat("\n====================\n")
cat("Scenario:", scn$name, "\n")
cat("====================\n\n")

cat("Observed PK:\n")
cat("Cmax =", obs_cmax, "umol/L\n")
cat("Tmax =", obs_tmax_h, "h\n")
cat("AUC_last =", obs_auc_last_h, "umol*h/L\n")
cat("lambda_z =", obs_lambda_z_per_h, "1/h\n")
cat("AUC_inf =", obs_auc_inf_h, "umol*h/L\n")
cat("CL =", obs_cl_L_per_min_per_kg, "L/min/kg\n")
cat("Vd =", obs_vd_L_per_kg, "L/kg\n")
cat("adj_r2 =", lz_obs$adj_r2, "\n")
cat("used terminal points =", paste(lz_obs$used_idx, collapse = ", "), "\n\n")

cat("Simulated PK from OSP Suite:\n")
cat("C_max =", sim_cmax, sim_cmax_unit, "\n")
cat("t_max =", sim_tmax_h, "h", "(raw:", sim_tmax, sim_tmax_unit, ")\n")
cat("AUC_last_obsWindow =", sim_auc_last_h, "umol*h/L", "(raw:", sim_auc_last, sim_auc_last_unit, ")\n")
cat("AUC_inf =", sim_auc_inf_h, "umol*h/L", "(raw:", sim_auc_inf, sim_auc_inf_unit, ")\n")
cat("CL =", sim_cl_L_per_min_per_kg, "L/min/kg", "(raw:", sim_cl, sim_cl_unit, ")\n")
cat("Vd =", sim_vd_L_per_kg, "L/kg", "(raw:", sim_vd, sim_vd_unit, ")\n\n")

cat("Comparison table:\n")
print(comparison)

cat("\nTmax comparison:\n")
print(tmax_table)

cat("\nOverall GMFE =", overall_gmfe, "\n")

# -----------------------------
# optional plot
# -----------------------------
dc <- DataCombined$new()
dc$addDataSets(obs_data)
dc$addSimulationResults(
  simulationResults = res,
  quantitiesOrPaths = out_path,
  names = scn$name
)

plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)