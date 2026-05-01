library(ospsuite)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)

simulation_files <- c(
  "~/Desktop/Thesis/TMDD/Mouse TvcMMAE 10 mg_kg with TMDD and Tumor.pkml"
  
)

output_paths <- c(
  "Organism|VenousBlood|Plasma|TvcMMAE|conjugatedMMAE",
  "Organism|VenousBlood|Plasma|nAb|totalmAb",
  "Organism|Tumor|TvcMMAE|conjugatedMMAEtumor",
  "Organism|Tumor|nAb|totalmAbtumor"
)

pk_parameters <- c(
  "AUC_tEnd",
  "AUC_inf",
  "C_max",
  "t_max",
  "Thalf"
)

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

output_file <- file.path(
  "~/Desktop/Thesis/Results",
  paste0("pk-summary-", timestamp, ".csv")
)


extract_pk_parameters <- function(simulation_file,
                                  output_paths,
                                  pk_parameters) {
  simulation <- loadSimulation(simulation_file)
  
  simulation_results <- runSimulations(
    simulations = simulation
  )[[1]]
  
  pk_analysis <- calculatePKAnalyses(
    results = simulation_results
  )
  
  pkAnalysesToDataFrame(pk_analysis) |>
    filter(
      QuantityPath %in% output_paths,
      Parameter %in% pk_parameters
    ) |>
    mutate(
      simulation_file = basename(simulation_file),
      Parameter = recode(
        Parameter,
        AUC_tEnd = "AUCtend",
        AUC_inf = "AUCtinf",
        C_max = "Cmax",
        t_max = "Tmax"
      )
    ) |>
    select(
      simulation_file,
      IndividualId,
      QuantityPath,
      Parameter,
      Value,
      Unit
    ) |>
    pivot_wider(
      names_from = Parameter,
      values_from = c(Value, Unit),
      names_glue = "{Parameter}_{.value}"
    )
}

pk_summary <- map_dfr(
  simulation_files,
  extract_pk_parameters,
  output_paths = output_paths,
  pk_parameters = pk_parameters
)

write_csv(pk_summary, output_file)

pk_summary
