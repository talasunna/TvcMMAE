library(ospsuite)
library(ggplot2)

options(width = 140)
options(ospsuite.plots.watermarkEnabled = FALSE)

# User inputs -------------------------------------------------------------

mobi_path <- "~/Desktop/Thesis/TMDD/Mouse TvcMMAE 10 mg_kg with TMDD and Tumor.pkml"
include_observed <- TRUE

observed_spec <- data.frame(
  dataset_path = c(
    #"~/Desktop/Thesis/Datasets/Observed_data/plasmaPK/Chang et al.FreeMMAE.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/plasmaPK/Chang et al.conjMMAE.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/plasmaPK/Chang et al.totalMMAE.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/plasmaPK/Chang et al.TvcMMAE.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/tumorPK/Chang et al.freeMMAE.tumor.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/tumorPK/Chang et al.conjMMAE.tumor.pkml"
    #"~/Desktop/Thesis/Datasets/Observed_data/tumorPK/Chang et al.totalMMAE.tumor.pkml"
    "~/Desktop/Thesis/Datasets/Observed_data/tumorPK/Chang et al.TvcMMAE.tumor.pkml"
    
    
  ),
  label = c(
    #"Observed Free MMAE in Tumor"
    #"Observed Conjugated MMAE in Tumor"
    #"Observed Total MMAE in Tumor"
    "Observed Total mAb in Tumor"
  ),
  color = c(
    "#0072B2"
    #"darkgreen"
    #"red"
  ),
  shape = c(
    16
    #17,
    #18
  ),
  stringsAsFactors = FALSE
)

plot_title <- "Total mAb Concentration in Tumor vs. Time"

curve_spec <- data.frame(
  sim_path = c(
    #"Organism|Tumor|MMAE|freeMMAEtumor"
    #"Organism|Tumor|TvcMMAE|conjugatedMMAEtumor"
    #"Organism|Tumor|MMAE|totalMMAEtumor"
    "Organism|Tumor|nAb|totalmAbtumor"
    
    
    
  ),
  label = c(
    #"Free MMAE Concentration in Tumor"
    #"Conjugated MMAE in Tumor"
    #"Total MMAE in Tumor"
    "Total mAb in Tumor"
    
  ),
  color = c(
    "#E69F00" 
    #"darkgreen"
    #"red"
  ),
  linewidth = c(
    0.8
    #0.8
    #0.8
  ),
  linetype = c(
    "solid"
    #"dashed"
    #"dashed"
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
  
  if (include_observed && length(observed_data) > 0) {
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
    obs_linetype_values <- stats::setNames(
      rep("blank", nrow(observed_spec)),
      observed_spec$label
    )
    obs_linewidth_values <- stats::setNames(
      rep(0, nrow(observed_spec)),
      observed_spec$label
    )
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
    title,
    include_observed = FALSE,
    x_unit = ospUnits$Time$h,
    y_scale = "log"
) {
  base_plot <- ospsuite::plotTimeProfile(
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
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.60, 0.40),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "grey70"),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

plot_pk_profile <- function(
    mobi_path,
    curve_spec,
    plot_title,
    include_observed = FALSE,
    observed_spec = NULL,
    x_unit = ospUnits$Time$h,
    y_scale = "log"
) {
  simulation_results <- load_simulation_results(mobi_path)
  
  observed_data <- NULL
  if (include_observed) {
    observed_data <- load_observed_data(observed_spec)
  }
  
  plot_data <- build_plot_data(
    simulation_results = simulation_results,
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
  observed_spec = observed_spec
)

print(p)