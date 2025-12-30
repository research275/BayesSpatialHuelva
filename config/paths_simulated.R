# =============================================================================
# PATH CONFIGURATION FOR SIMULATED DATA
# =============================================================================
# This file configures paths to use simulated data instead of real data.
# To use simulated data, source this file instead of config/paths.R
# =============================================================================

BASE_DIR <- here::here()
DATA_DIR <- file.path(BASE_DIR, "data", "simulated")
RESULTS_DIR <- file.path(BASE_DIR, "results")
SRC_DIR <- file.path(BASE_DIR, "src")

PATHS <- list(
  boundaries = file.path(DATA_DIR, "boundaries", "Malla_sim", "malla.gpkg"),
  
  samples = file.path(DATA_DIR, "samples", "Datos", "HuelvaSimulado.gpkg"),
  
  prediction_grid = file.path(DATA_DIR, "prediction_grid", "Grid_sim", "capa_sin_agua.gpkg"),
  
  land_use_samples = list(
    agric = file.path(DATA_DIR, "land_use", "samples", "Agric", "Agric.gpkg"),
    urban = file.path(DATA_DIR, "land_use", "samples", "Urban", "Urban.gpkg"),
    industria = file.path(DATA_DIR, "land_use", "samples", "Industria", "Industria.gpkg"),
    refineria = file.path(DATA_DIR, "land_use", "samples", "Refineria", "Refineria.gpkg"),
    phospho = file.path(DATA_DIR, "land_use", "samples", "Phospho", "Phospho.gpkg"),
    marsh = file.path(DATA_DIR, "land_use", "samples", "Marsh", "Marsh.gpkg"),
    bare = file.path(DATA_DIR, "land_use", "samples", "Bare", "Bare.gpkg"),
    park = file.path(DATA_DIR, "land_use", "samples", "Park", "Park.gpkg")
  ),
  
  land_use_covariates = list(
    agric = file.path(DATA_DIR, "land_use", "covariates", "Agric", "Agric.gpkg"),
    urban = file.path(DATA_DIR, "land_use", "covariates", "Urban", "Urban.gpkg"),
    industria = file.path(DATA_DIR, "land_use", "covariates", "Industria", "Industria.gpkg"),
    refineria = file.path(DATA_DIR, "land_use", "covariates", "Refineria", "Refineria.gpkg"),
    phospho = file.path(DATA_DIR, "land_use", "covariates", "Phospho", "Phospho.gpkg"),
    marsh = file.path(DATA_DIR, "land_use", "covariates", "Marsh", "Marsh.gpkg"),
    bare = file.path(DATA_DIR, "land_use", "covariates", "Bare", "Bare.gpkg"),
    park = file.path(DATA_DIR, "land_use", "covariates", "Park", "Park.gpkg")
  ),
  
  land_use_polygons = list(
    agric = file.path(DATA_DIR, "land_use", "polygons", "Agric.gpkg"),
    urban = file.path(DATA_DIR, "land_use", "polygons", "Urban.gpkg"),
    industria = file.path(DATA_DIR, "land_use", "polygons", "Industria.gpkg"),
    refineria = file.path(DATA_DIR, "land_use", "polygons", "Refineria.gpkg"),
    phospho = file.path(DATA_DIR, "land_use", "polygons", "Phospho.gpkg"),
    marsh = file.path(DATA_DIR, "land_use", "polygons", "Marsh.gpkg"),
    bare = file.path(DATA_DIR, "land_use", "polygons", "Bare.gpkg"),
    park = file.path(DATA_DIR, "land_use", "polygons", "Park.gpkg")
  )
)

OUTPUT_PATHS <- list(
  figures = file.path(RESULTS_DIR, "figures"),
  reports = file.path(RESULTS_DIR, "reports"),
  models = file.path(RESULTS_DIR, "models"),
  predictions = file.path(RESULTS_DIR, "predictions")
)

create_output_dirs <- function() {
  lapply(OUTPUT_PATHS, function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      cat("Created directory:", path, "\n")
    }
  })
  invisible(TRUE)
}

check_input_files <- function() {
  missing_files <- c()
  
  individual_files <- list(
    PATHS$boundaries,
    PATHS$samples,
    PATHS$prediction_grid
  )
  
  for (file in individual_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  for (file in PATHS$land_use_samples) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  for (file in PATHS$land_use_covariates) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("WARNING: The following files were not found:\n")
    for (file in missing_files) {
      cat("  -", file, "\n")
    }
    return(FALSE)
  } else {
    cat("All input files (simulated) are available.\n")
    return(TRUE)
  }
}

cat("Path configuration for SIMULATED DATA loaded.\n")
cat("Use create_output_dirs() to create output directories.\n")
cat("Use check_input_files() to verify input files.\n")

