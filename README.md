# BayesSpatialHuelva

## What is this repository for?

BayesSpatialHuelva is a software for Bayesian spatial modeling and prediction of heavy metal contamination (Copper and Lead) in soil samples. It implements and compares three inference approaches: INLA (Integrated Nested Laplace Approximations), MCMC (Markov Chain Monte Carlo), and Variational Bayes (ADVI). The software computes posterior distributions for geospatial regression models with Matern covariance structures and generates exceedance probability maps.

This repository contains simulated data similar to the one used in the manuscript due to restrictions on the original dataset.

## How do I get set up?

### Requirements

- R version 4.0 or higher
- The following R packages:
  - `sf`, `dplyr`, `tidyr`, `Matrix`, `fields`
  - `INLA` (install from r-inla.org)
  - `spBayes`, `rstan`
  - `ggplot2`, `viridis`, `cowplot`, `plotly`
  - `here`, `writexl`, `readxl`

### Installation

1. Clone or download this repository
2. Install the required R packages
3. For INLA, run: `install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)`
4. Configure the paths in `config/paths.R` or use `config/paths_simulated.R` for the simulated dataset

## Usage

### Main Scripts

1. **algorithms_cu.R** - Bayesian spatial analysis for Copper (Cu) contamination

   - Loads sample data and land use covariates
   - Fits spatial models using INLA, MCMC, and Variational Bayes
   - Generates predictions on a spatial grid
2. **algorithms_pb.R** - Bayesian spatial analysis for Lead (Pb) contamination

   - Same methodology as Cu analysis but for Lead concentrations
3. **exceedance_plots.R** - Generates exceedance probability maps

   - Compares exceedance probabilities across the three inference methods
   - Produces visualization outputs
4. **computation.R** - Simulation study script

   - Runs Monte Carlo, INLA and ADVI experiments under different scenarios
   - Compares computational performance and accuracy of methods

### Running the Analysis

1. Configure the data paths in `config/paths.R`
2. Run the desired analysis script in R:
   ```r
   source("algorithms_cu.R")  # For Copper analysis
   source("algorithms_pb.R")  # For Lead analysis
   source("exceedance_plots.R")  # For exceedance probability maps
   ```

### Data Structure

The repository expects the following data structure:

- `data/simulated/boundaries/` - Study area boundary files
- `data/simulated/samples/` - Sample point locations with metal concentrations
- `data/simulated/prediction_grid/` - Grid for spatial predictions
- `data/simulated/land_use/` - Land use polygons and covariates (Agricultural, Urban, Industrial, Refinery, Phosphogypsum, Marsh, Bare, Park)

### Output Files

Results are saved in the `results/` directory:

- `results/figures/` - Generated plots and maps
- `results/models/` - Fitted model objects
- `results/predictions/` - Prediction grids and exceedance probabilities
- `results/reports/` - Summary statistics and comparison tables

## License

Copyright (C) 2026 Anonymous Authors. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
