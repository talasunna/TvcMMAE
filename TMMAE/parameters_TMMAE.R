library(ospsuite)
library(ggplot2)

rm(list = ls())

# --- paths ---
sim_path <- "~/Desktop/Thesis/TMMAE/Mouse1 _DC_i.v. _10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/TMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# optional: inspect what was loaded
print(obs_data)

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
  list(name = "FcRnKd: 1.07 μmol/l, HER2Kd: 3.92 μmol/l, HER2 ref: 1 μmol/l, Dose: 15 mg/kg, FU = 0.5", fcrn_endo_kd = 1, tmmae_her2_kd = 3.92, herref = 1, dose = NULL, fu = 0.2)
  #list(name = "Dose: 15 mg/kg" ,             fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.5e-5),
  #list(name = "Dose: 16 mg/kg",             fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.6e-5),
  #list(name = "Dose: 17 mg/kg",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.7e-5),
  #list(name = "Dose: 18 mg/kg",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.8e-5),
  #list(name = "Dose: 19 mg/kg",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.9e-5)
  #list(name = "fu = 1.5",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL, fu = 1.5),
  #list(name = "fu = 0.8",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL, fu = 0.8), 
  #list(name = "fu = 0.5",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL, fu = 0.5),
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
}

# --- plot in hours ---
plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(dc, plot_cfg)