library(ospsuite)
library(ggplot2)

options(width = 140)
options(ospsuite.plots.watermarkEnabled = FALSE)

# User inputs -------------------------------------------------------------

mobi_path <- "~/Desktop/Thesis/TMDD/Mouse TvcMMAE 10 mg_kg with TMDD and Tumor.pkml"

include_observed <- TRUE
dataset_path <- "~/Desktop/Thesis/MMAE/Chang et al.MMAE.IV__10mgKg_MMAE_Mouse.pkml"

plot_title <- "MMAE concentration-time profiles"
obs_label <- "Observed data"
obs_color <- "orange"
obs_shape <- 16

curve_spec <- data.frame(
  sim_path = c(
    "Organism|VenousBlood|Plasma|MMAE|Concentration in container",
    "Organism|Tumor|Interstitial|MMAE|Concentration in container",
    "Organism|Tumor|Intracellular|MMAE|Concentration in container",
    "Organism|Tumor|Interstitial|nAb|Concentration in container"
    
  ),
  label = c(
    "Plasma MMAE",
    "Tumor interstitial MMAE",
    "Tumor intracellular MMAE",
    "Tumor interstitial nAb"
    
  ),
  color = c(
    "red",
    "darkgreen",
    "blue",
    "pink"
  ),
  linewidth = c(
    0.8,
    0.8,
    0.8,
    0.8
  ),
  linetype = c(
    "solid",
    "solid",
    "solid",
    "solid"
    
  ),
  stringsAsFactors = FALSE
)

# Helper functions --------------------------------------------------------

load_simulation_results <- function(file_path) {
  simulation <- loadSimulation(file_path)
  runSimulations(simulation)[[1]]
}

load_observed_data <- function(file_path, label) {
  data_set <- loadDataSetFromPKML(filePath = file_path)
  data_set$name <- label
  data_set
}

build_plot_data <- function(
    simulation_results,
    curve_spec,
    include_observed = FALSE,
    observed_data = NULL
) {
  plot_data <- DataCombined$new()
  
  for (i in seq_len(nrow(curve_spec))) {
    plot_data$addSimulationResults(
      simulationResults = simulation_results,
      quantitiesOrPaths = curve_spec$sim_path[i],
      names = curve_spec$label[i]
    )
  }
  
  if (include_observed) {
    plot_data$addDataSets(observed_data)
  }
  
  plot_data
}

build_legend_spec <- function(
    curve_spec,
    include_observed = FALSE,
    obs_label = "Observed data",
    obs_color = "red",
    obs_shape = 16
) {
  color_values <- stats::setNames(curve_spec$color, curve_spec$label)
  linetype_values <- stats::setNames(curve_spec$linetype, curve_spec$label)
  linewidth_values <- stats::setNames(curve_spec$linewidth, curve_spec$label)
  shape_values <- stats::setNames(rep(NA, nrow(curve_spec)), curve_spec$label)
  
  if (include_observed) {
    color_values <- c(
      color_values,
      stats::setNames(obs_color, obs_label)
    )
    linetype_values <- c(
      linetype_values,
      stats::setNames("blank", obs_label)
    )
    linewidth_values <- c(
      linewidth_values,
      stats::setNames(0, obs_label)
    )
    shape_values <- c(
      shape_values,
      stats::setNames(obs_shape, obs_label)
    )
    legend_order <- c(curve_spec$label, obs_label)
  } else {
    legend_order <- curve_spec$label
  }
  
  list(
    order = legend_order,
    color_values = color_values,
    linetype_values = linetype_values,
    linewidth_values = linewidth_values,
    shape_values = shape_values
  )
}

make_time_profile_plot <- function(
    plot_data,
    legend_spec,
    title,
    include_observed = FALSE,
    x_unit = ospUnits$Time$h,
    y_scale = "log"
) {
  base_plot <- plotTimeProfile(
    plotData = plot_data,
    xUnit = x_unit,
    yScale = y_scale,
    mapping = aes(
      color = name,
      linetype = name,
      linewidth = name
    ),
    observedMapping = if (include_observed) {
      aes(color = name, shape = name)
    } else {
      NULL
    },
    showLegendPerDataset = "all"
  )
  
  base_plot +
    scale_color_manual(
      values = legend_spec$color_values,
      breaks = legend_spec$order
    ) +
    scale_linetype_manual(
      values = legend_spec$linetype_values,
      breaks = legend_spec$order
    ) +
    scale_linewidth_manual(
      values = legend_spec$linewidth_values,
      breaks = legend_spec$order
    ) +
    scale_shape_manual(
      values = legend_spec$shape_values,
      breaks = legend_spec$order
    ) +
    labs(
      title = title,
      x = "Time [h]",
      y = "Concentration [\u03bcmol/l]",
      color = NULL,
      linetype = NULL,
      linewidth = NULL,
      shape = NULL
    ) +
    guides(
      linetype = "none",
      linewidth = "none",
      shape = "none",
      color = guide_legend(
        override.aes = list(
          linetype = unname(legend_spec$linetype_values[legend_spec$order]),
          linewidth = unname(legend_spec$linewidth_values[legend_spec$order]),
          shape = unname(legend_spec$shape_values[legend_spec$order])
        )
      )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank()
    )
}

plot_pk_profile <- function(
    mobi_path,
    curve_spec,
    plot_title,
    include_observed = FALSE,
    dataset_path = NULL,
    obs_label = "Observed data",
    obs_color = "red",
    obs_shape = 16,
    x_unit = ospUnits$Time$h,
    y_scale = "log"
) {
  simulation_results <- load_simulation_results(mobi_path)
  
  observed_data <- NULL
  if (include_observed) {
    observed_data <- load_observed_data(dataset_path, obs_label)
  }
  
  plot_data <- build_plot_data(
    simulation_results = simulation_results,
    curve_spec = curve_spec,
    include_observed = include_observed,
    observed_data = observed_data
  )
  
  legend_spec <- build_legend_spec(
    curve_spec = curve_spec,
    include_observed = include_observed,
    obs_label = obs_label,
    obs_color = obs_color,
    obs_shape = obs_shape
  )
  
  make_time_profile_plot(
    plot_data = plot_data,
    legend_spec = legend_spec,
    title = plot_title,
    include_observed = include_observed,
    x_unit = x_unit,
    y_scale = y_scale
  )
}

# Run ---------------------------------------------------------------------

p <- plot_pk_profile(
  mobi_path = mobi_path,
  curve_spec = curve_spec,
  plot_title = plot_title,
  include_observed = include_observed,
  dataset_path = dataset_path,
  obs_label = obs_label,
  obs_color = obs_color,
  obs_shape = obs_shape
)

print(p)