library(ospsuite)
library(ospsuite.reportingengine)

rm(list = ls())

sim_path <- "~/Desktop/Thesis/Mouse MMAE 0.10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|MMAE|Concentration in container"
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.FreeMMAE.pkml"

auc_end_time_h <- 168
auc_end_time_min <- auc_end_time_h * 60

calc_auc_to_time <- function(time_min, conc, end_time_min) {
  keep <- is.finite(time_min) & is.finite(conc) & time_min <= end_time_min
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

obs_data <- loadDataSetFromPKML(filePath = dataset_path)
obs_df <- dataSetToDataFrame(obs_data)

obs_auc <- calc_auc_to_time(
  time_min = obs_df$xValues * 60,
  conc = obs_df$yValues,
  end_time_min = auc_end_time_min
)

sim <- loadSimulation(sim_path)
setOutputs(out_path, sim)

sim_results <- runSimulations(simulations = sim)[[1]]

sim_values <- getOutputValues(
  simulationResults = sim_results,
  quantitiesOrPaths = out_path
)

sim_df <- sim_values$data

sim_auc <- calc_auc_to_time(
  time_min = sim_df$Time,
  conc = sim_df[[out_path]],
  end_time_min = auc_end_time_min
)

gmfe <- ospsuite.reportingengine::calculateGMFE(
  x = sim_auc,
  y = obs_auc
)

error_summary <- data.frame(
  auc_end_time_h = auc_end_time_h,
  observed_auc_umol_min_l = round(obs_auc, 3),
  predicted_auc_umol_min_l = round(sim_auc, 3),
  gmfe = round(gmfe, 3)
)

print(error_summary, row.names = FALSE)