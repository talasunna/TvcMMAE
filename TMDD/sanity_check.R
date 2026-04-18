library(ospsuite)

rm(list = ls())
options(width = 140)

# -------------------------
# file paths
# -------------------------
pksim_path   <- "~/Desktop/Thesis/Mouse TvcMMAE i.v. 10 mg_kg.pkml"
mobi_path    <- "~/Desktop/Thesis/TMDD/Mouse TvcMMAE i.v. 10 mg_kg with TMDD and Tumor.pkml"
dataset_path <- "~/Desktop/Thesis/TMMAE/Chang_et_al.PK_observed.IV__10mgKg_ADC_Mouse.pkml"

# output path
out_path <- "Organism|VenousBlood|Plasma|TvcMMAE|Concentration in container"

# -------------------------
# load simulations
# -------------------------
sim_pksim <- loadSimulation(pksim_path)
sim_mobi  <- loadSimulation(mobi_path)

# -------------------------
# run simulations
# -------------------------
res_pksim <- runSimulations(sim_pksim)[[1]]
res_mobi  <- runSimulations(sim_mobi)[[1]]

# -------------------------
# load observed data
# -------------------------
obs_data <- loadDataSetFromPKML(filePath = dataset_path)
obs_data$name <- "Observed data"

# -------------------------
# combine all data
# -------------------------
dc <- DataCombined$new()

dc$addSimulationResults(
  simulationResults = res_pksim,
  quantitiesOrPaths = out_path,
  names = "PK"
)

dc$addSimulationResults(
  simulationResults = res_mobi,
  quantitiesOrPaths = out_path,
  names = "PK with TMDD and Tumor"
)

dc$addDataSets(obs_data)

# -------------------------
# plot
# -------------------------
plot_cfg <- DefaultPlotConfiguration$new()
plot_cfg$xUnit <- ospUnits$Time$h
plot_cfg$yAxisScale <- tlf::Scaling$log

plotIndividualTimeProfile(
  dataCombined = dc,
  defaultPlotConfiguration = plot_cfg
)