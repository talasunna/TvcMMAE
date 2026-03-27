library(ospsuite)

rm(list = ls())

# --- paths ---
sim_path <- "~/Desktop/Thesis/Mouse1 _DC_i.v. _10 mg_kg.pkml"
out_path <- "Organism|VenousBlood|Plasma|TMMAE|Concentration"
dataset_path <- "~/Desktop/Thesis/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# --- load observed data ---
obs_data <- loadDataSetFromPKML(filePath = dataset_path)

# optional: inspect what was loaded
print(obs_data)

# --- function that applies one scenario to a simulation ---
apply_changes <- function(sim,
                          fcrn_endo_kd = NULL,
                          tmmae_her2_kd = NULL,
                          herref = NULL,
                          dose = NULL) {
  
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
  #list(name = "Baseline simulation", fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = NULL),
  list(name = "dose 10",             fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 1.00e-05),
  list(name = "dose 20",             fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 2.00e-05)
  #list(name = "dose 5",              fcrn_endo_kd = NULL, tmmae_her2_kd = NULL, herref = NULL, dose = 5)
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
    dose          = scn$dose
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