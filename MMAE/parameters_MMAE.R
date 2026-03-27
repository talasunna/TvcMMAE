library(ospsuite)
library(ggplot2)

rm(list = ls())

# --- paths ---
sim_path <- "~/Desktop/Thesis/MMAE/Mouse1 MMAE i.v. 10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|MMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.MMAE.IV__10mgKg_MMAE_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# optional: inspect what was loaded
print(obs_data)

# --- function that applies one scenario to a simulation ---
apply_changes <- function(sim,
                          fu = NULL,
                          lipophilicity = NULL,
                          specific_clearance = NULL) {
  
  if (!is.null(fu)) {
    fu_plasma <- getParameter("MMAE|Fraction unbound (plasma, reference value)", sim)
    setParameterValues(fu_plasma, fu)
    
    if (!is.null(lipophilicity)) {
      lip <- getParameter("MMAE|Lipophilicity", sim)
      setParameterValues(lip, lipophilicity)
    }
    
    if (!is.null(specific_clearance)) {
      sc<- getParameter("MMAE-Total Hepatic Clearance-Mouse|Specific clearance", sim)
      setParameterValues(sc, specific_clearance)
    }
  }

  
  sim$solver$absTol <- 1e-12
  sim$solver$relTol <- 1e-8
  
  invisible(sim)
}

# --- scenarios ---
scenarios <- list(
  list(name = "FU = 0.5, Lip = 3, sc = 1.31", specific_clearance = 0.02, lipophilicity = 2.7, fu = 0.6)
)

# --- combine observed + simulated data ---
dc <- DataCombined$new()

# add observed data 
dc$addDataSets(obs_data)

for (scn in scenarios) {
  sim <- loadSimulation(sim_path)
  setOutputs(out_path, sim)
  
  apply_changes(
    sim,
    fu            = scn$fu,
    lipophilicity = scn$lipophilicity,
    specific_clearance = scn$specific_clearance
  )
  
  res <- runSimulations(sim)[[1]]
  
  dc$addSimulationResults(
    simulationResults = res,
    quantitiesOrPaths = out_path,
    names = scn$name
  )
}

# --- plot in hours ---
plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)