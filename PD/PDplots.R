library(ospsuite)
library(ggplot2)

options(width = 140)
options(ospsuite.plots.watermarkEnabled = FALSE)

# User inputs -------------------------------------------------------------

include_observed <- TRUE

# Observed data: put control here, and optionally the observed dose data too
observed_spec <- data.frame(
  dataset_path = c(
    "~/Desktop/Thesis/Data/Observed_data/PD/Chang et al.TGI.ctrl.pkml"
    #"~/Desktop/Thesis/Data/Observed_data/PD/Chang et al.TGI.1mgkg.pkml"
    # "~/Desktop/Thesis/Data/Observed_data/PD/Chang et al.TGI.3mgkg.pkml",
    # "~/Desktop/Thesis/Data/Observed_data/PD/Chang et al.TGI.10mgkg.pkml"
  ),
  label = c(
    "Observed, control"
    #"Observed, 1 mg/kg"
    # "Observed, 3 mg/kg",
    # "Observed, 10 mg/kg"
  ),
  color = c(
    #"#000000",
    "#0072B2"
    # "#009E73",
    # "#D55E00"
  ),
  shape = c(
    16
    #16
    # 16,
    # 16
  ),
  stringsAsFactors = FALSE
)

# Simulation curves: one row per simulation
curve_spec <- data.frame(
  mobi_path = c(
    "~/Desktop/Thesis/PD/TGI 1 mg_kg.pkml",
    "~/Desktop/Thesis/PD/TGI 3 mg_kg.pkml",
    "~/Desktop/Thesis/PD/TGI 10 mg_kg-2.pkml"
  ),
  sim_path = c(
    "Organism|Tumor|Volume",
    "Organism|Tumor|Volume",
    "Organism|Tumor|Volume"
  ),
  label = c(
    "Predicted, 1 mg/kg",
    "Predicted, 3 mg/kg",
    "Predicted, 10 mg/kg"
  ),
  color = c(
    "#CC79A7",
    "#009E73",
    "#D55E00"
  ),
  linewidth = c(
    0.8,
    0.8,
    0.8
  ),
  linetype = c(
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

load_observed_data <- function(observed_spec) {
  observed_data <- vector("list", nrow(observed_spec))
  
  for (i in seq_len(nrow(observed_spec))) {
    data_set <- loadDataSetFromPKML(filePath = observed_spec$dataset_path[i])
    data_set$name <- observed_spec$label[i]
    observed_data[[i]] <- data_set
  }
  
  observed_data
}

build_plot_data <- function(
    curve_spec,
    include_observed = FALSE,
    observed_data = NULL
) {
  plot_data <- DataCombined$new()
  
  for (i in seq_len(nrow(curve_spec))) {
    simulation_results <- load_simulation_results(curve_spec$mobi_path[i])
    
    plot_data$addSimulationResults(
      simulationResults = simulation_results,
      quantitiesOrPaths = curve_spec$sim_path[i],
      names = curve_spec$label[i]
    )
  }
  
  if (include_observed && !is.null(observed_data) && length(observed_data) > 0) {
    plot_data$addDataSets(observed_data)
  }
  
  plot_data
}

build_legend_spec <- function(
    curve_spec,
    observed_spec = NULL,
    include_observed = FALSE
) {
  color_values <- stats::setNames(curve_spec$color, curve_spec$label)
  linetype_values <- stats::setNames(curve_spec$linetype, curve_spec$label)
  linewidth_values <- stats::setNames(curve_spec$linewidth, curve_spec$label)
  shape_values <- stats::setNames(rep(NA, nrow(curve_spec)), curve_spec$label)
  
  legend_order <- curve_spec$label
  
  if (include_observed && !is.null(observed_spec) && nrow(observed_spec) > 0) {
    obs_color_values <- stats::setNames(observed_spec$color, observed_spec$label)
    obs_linetype_values <- stats::setNames(rep("blank", nrow(observed_spec)), observed_spec$label)
    obs_linewidth_values <- stats::setNames(rep(0, nrow(observed_spec)), observed_spec$label)
    obs_shape_values <- stats::setNames(observed_spec$shape, observed_spec$label)
    
    color_values <- c(color_values, obs_color_values)
    linetype_values <- c(linetype_values, obs_linetype_values)
    linewidth_values <- c(linewidth_values, obs_linewidth_values)
    shape_values <- c(shape_values, obs_shape_values)
    
    legend_order <- c(curve_spec$label, observed_spec$label)
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
    include_observed = FALSE,
    x_unit = ospUnits$Time$d,
    y_unit = "µl",
    y_scale = "linear"
) {
  base_plot <- ospsuite::plotTimeProfile(
    plotData = plot_data,
    xUnit = x_unit,
    yUnit = y_unit,
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
    showLegendPerDataset = "observed"
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
      x = "Time [days]",
      y = expression("Tumor volume ["*mu*"l]"),
      color = NULL,
      linetype = NULL,
      linewidth = NULL,
      shape = NULL
    ) +
    guides(
      linetype = "none",
      linewidth = "none",
      shape = "none",
      color = guide_legend(ncol = 1)
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 9),
      plot.margin = margin(5.5, 15, 5.5, 5.5)
    )
}

plot_tgi_profiles <- function(
    curve_spec,
    include_observed = FALSE,
    observed_spec = NULL,
    x_unit = ospUnits$Time$d,
    y_unit = "µl",
    y_scale = "linear"
) {
  observed_data <- NULL
  
  if (include_observed) {
    observed_data <- load_observed_data(observed_spec)
  }
  
  plot_data <- build_plot_data(
    curve_spec = curve_spec,
    include_observed = include_observed,
    observed_data = observed_data
  )
  
  legend_spec <- build_legend_spec(
    curve_spec = curve_spec,
    observed_spec = observed_spec,
    include_observed = include_observed
  )
  
  make_time_profile_plot(
    plot_data = plot_data,
    legend_spec = legend_spec,
    include_observed = include_observed,
    x_unit = x_unit,
    y_unit = y_unit,
    y_scale = y_scale
  )
}

# Run ---------------------------------------------------------------------

p <- plot_tgi_profiles(
  curve_spec = curve_spec,
  include_observed = include_observed,
  observed_spec = observed_spec,
  x_unit = ospUnits$Time$d,
  y_unit = "µl"
)

print(p)