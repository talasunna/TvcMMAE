library(ospsuite)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

simulation_files <- c(
  "~/Desktop/Thesis/TMDD/Mouse TvcMMAE 10 mg_kg with TMDD and Tumor.pkml"
)

output_paths <- c(
  "Organism|VenousBlood|Plasma|TvcMMAE|conjugatedMMAE",
  "Organism|VenousBlood|Plasma|MMAE|totalMMAEplasma",
  "Organism|VenousBlood|Plasma|nAb|totalmAbplasma",
  "Organism|Tumor|MMAE|freeMMAEtumor",
  "Organism|Tumor|TvcMMAE|conjugatedMMAEtumor",
  "Organism|Tumor|MMAE|totalMMAEtumor",
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

output_dir <- "~/Desktop/Thesis/Results"

pk_output_file <- file.path(
  output_dir,
  paste0("pk-summary-", timestamp, ".csv")
)

gmfe_output_file <- file.path(
  output_dir,
  paste0("pk-gmfe-table-", timestamp, ".csv")
)

gmfe_summary_file <- file.path(
  output_dir,
  paste0("pk-gmfe-summary-", timestamp, ".csv")
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

write_csv(pk_summary, pk_output_file)


# =============================================================================
# GMFE CALCULATION
# =============================================================================
# Observed IV values from Chang et al. Table 2
#
# Units:
# Tmax    = h
# Cmax    = nM
# AUCtend = h*nM
# AUCtinf = h*nM
# =============================================================================

observed_iv_pk <- tribble(
  ~Analyte,             ~Matrix,  ~Tmax,  ~Cmax, ~AUCtend, ~AUCtinf,
  
  "Total antibody",     "Plasma", 0.167,  1690,  1.34e5,   7.48e5,
  "Total antibody",     "Tumor",  72,     242,   3.52e4,   1.94e5,
  
  "Total MMAE",         "Plasma", 0.167,  4058,  1.28e5,   1.52e5,
  "Total MMAE",         "Tumor",  72,     693,   8.79e4,   1.74e5,
  
  "Unconjugated MMAE",  "Plasma", 0.167,  11.9,  224,      324,
  "Unconjugated MMAE",  "Tumor",  72,     419,   4.79e4,   8.16e4,
  
  "Conjugated MMAE",    "Plasma", 0.167,  4046,  1.28e5,   1.51e5,
  "Conjugated MMAE",    "Tumor",  24,     353,   4.01e4,   5.25e4
) |>
  pivot_longer(
    cols = c(Tmax, Cmax, AUCtend, AUCtinf),
    names_to = "PKParameter",
    values_to = "Observed"
  )


# Map your simulated output paths to the observed analytes from the paper

path_map <- tribble(
  ~QuantityPath,                                                    ~Analyte,            ~Matrix,
  
  "Organism|VenousBlood|Plasma|nAb|totalmAbplasma",                 "Total antibody",    "Plasma",
  "Organism|Tumor|nAb|totalmAbtumor",                               "Total antibody",    "Tumor",
  
  "Organism|VenousBlood|Plasma|MMAE|totalMMAEplasma",               "Total MMAE",        "Plasma",
  "Organism|Tumor|MMAE|totalMMAEtumor",                             "Total MMAE",        "Tumor",
  
  "Organism|VenousBlood|Plasma|MMAE|Concentration in container",    "Unconjugated MMAE", "Plasma",
  "Organism|Tumor|MMAE|freeMMAEtumor",                              "Unconjugated MMAE", "Tumor",
  
  "Organism|VenousBlood|Plasma|TvcMMAE|conjugatedMMAE",             "Conjugated MMAE",   "Plasma",
  "Organism|Tumor|TvcMMAE|conjugatedMMAEtumor",                     "Conjugated MMAE",   "Tumor"
)


# Unit conversion helpers

normalize_unit <- function(unit) {
  unit |>
    tolower() |>
    str_replace_all("μ", "µ")
}

convert_cmax_to_nM <- function(value, unit) {
  unit <- normalize_unit(unit)
  
  case_when(
    str_detect(unit, "µmol|umol") ~ value * 1000,
    str_detect(unit, "nmol|nm")   ~ value,
    TRUE ~ NA_real_
  )
}

convert_time_to_h <- function(value, unit) {
  unit <- normalize_unit(unit)
  
  case_when(
    str_detect(unit, "min") ~ value / 60,
    str_detect(unit, "h")   ~ value,
    TRUE ~ NA_real_
  )
}

convert_auc_to_h_nM <- function(value, unit) {
  unit <- normalize_unit(unit)
  
  concentration_factor <- case_when(
    str_detect(unit, "µmol|umol") ~ 1000,
    str_detect(unit, "nmol|nm")   ~ 1,
    TRUE ~ NA_real_
  )
  
  time_factor <- case_when(
    str_detect(unit, "min") ~ 1 / 60,
    str_detect(unit, "h")   ~ 1,
    TRUE ~ NA_real_
  )
  
  value * concentration_factor * time_factor
}


# Convert simulated PK values to the same units as the paper:
# Tmax = h, Cmax = nM, AUC = h*nM

simulated_pk_long <- pk_summary |>
  left_join(path_map, by = "QuantityPath") |>
  filter(!is.na(Analyte)) |>
  mutate(
    Tmax = convert_time_to_h(Tmax_Value, Tmax_Unit),
    Cmax = convert_cmax_to_nM(Cmax_Value, Cmax_Unit),
    AUCtend = convert_auc_to_h_nM(AUCtend_Value, AUCtend_Unit),
    AUCtinf = convert_auc_to_h_nM(AUCtinf_Value, AUCtinf_Unit)
  ) |>
  select(
    simulation_file,
    IndividualId,
    QuantityPath,
    Analyte,
    Matrix,
    Tmax,
    Cmax,
    AUCtend,
    AUCtinf
  ) |>
  pivot_longer(
    cols = c(Tmax, Cmax, AUCtend, AUCtinf),
    names_to = "PKParameter",
    values_to = "Predicted"
  )


# Calculate fold error and GMFE

gmfe_table <- simulated_pk_long |>
  inner_join(
    observed_iv_pk,
    by = c("Analyte", "Matrix", "PKParameter")
  ) |>
  filter(
    !is.na(Observed),
    !is.na(Predicted),
    Observed > 0,
    Predicted > 0
  ) |>
  mutate(
    Ratio_Predicted_Observed = Predicted / Observed,
    Fold_Error = pmax(Predicted / Observed, Observed / Predicted)
  ) |>
  arrange(
    Analyte,
    Matrix,
    PKParameter
  )


gmfe_summary <- gmfe_table |>
  summarise(
    GMFE = exp(mean(abs(log(Predicted / Observed)))),
    n_parameters = n()
  )


gmfe_by_analyte <- gmfe_table |>
  group_by(Analyte) |>
  summarise(
    GMFE = exp(mean(abs(log(Predicted / Observed)))),
    n_parameters = n(),
    .groups = "drop"
  )


gmfe_by_matrix <- gmfe_table |>
  group_by(Matrix) |>
  summarise(
    GMFE = exp(mean(abs(log(Predicted / Observed)))),
    n_parameters = n(),
    .groups = "drop"
  )


gmfe_by_parameter <- gmfe_table |>
  group_by(PKParameter) |>
  summarise(
    GMFE = exp(mean(abs(log(Predicted / Observed)))),
    n_parameters = n(),
    .groups = "drop"
  )


# Save results

write_csv(gmfe_table, gmfe_output_file)
write_csv(gmfe_summary, gmfe_summary_file)


# Print results in console

pk_summary

gmfe_table

gmfe_summary

gmfe_by_analyte

gmfe_by_matrix

gmfe_by_parameter