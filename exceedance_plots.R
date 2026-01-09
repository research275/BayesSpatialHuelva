# =============================================================================
# SPATIAL EXCEEDANCE PROBABILITY COMPARISON: INLA vs MCMC vs VB
# Lead (Pb) contamination in Huelva, Spain
# =============================================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(Matrix)
  library(fields)
  library(INLA)
  library(spBayes)
  library(rstan)
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(scales)
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# =============================================================================
# 1. LOAD DATA AND CONFIGURATION
# =============================================================================

cat("==========================================================\n")
cat("LOADING DATA\n")
cat("==========================================================\n")

source(here::here("config", "paths.R"))
check_input_files()
create_output_dirs()

crs_utm <- 25829

# Helper functions
read_sf_utm <- function(path, epsg) {
  x <- st_read(path, quiet = TRUE)
  if (is.na(st_crs(x))) stop(paste0("Missing CRS in ", path))
  if (st_crs(x)$epsg != epsg) x <- st_transform(x, epsg)
  x
}

min_dist_to <- function(points_sf, polys_sf) {
  as.numeric(st_distance(points_sf, st_union(polys_sf)))
}

# Load spatial data
boundary    <- read_sf_utm(PATHS$boundaries, crs_utm)
samples     <- read_sf_utm(PATHS$samples, crs_utm) |> st_zm()
pred_grid   <- read_sf_utm(PATHS$prediction_grid, crs_utm)

# Load land use polygons for covariates
land_use_samples <- list(
  Agric     = PATHS$land_use_samples$agric,
  Urban     = PATHS$land_use_samples$urban,
  Industry  = PATHS$land_use_samples$industria,
  Refinery  = PATHS$land_use_samples$refineria,
  Phospho   = PATHS$land_use_samples$phospho,
  Marsh     = PATHS$land_use_samples$marsh,
  Bare      = PATHS$land_use_samples$bare,
  Park      = PATHS$land_use_samples$park
)

land_use_grid <- list(
  Agric     = PATHS$land_use_covariates$agric,
  Urban     = PATHS$land_use_covariates$urban,
  Industry  = PATHS$land_use_covariates$industria,
  Refinery  = PATHS$land_use_covariates$refineria,
  Phospho   = PATHS$land_use_covariates$phospho,
  Marsh     = PATHS$land_use_covariates$marsh,
  Bare      = PATHS$land_use_covariates$bare,
  Park      = PATHS$land_use_covariates$park
)

land_use_polys_samples <- lapply(land_use_samples, read_sf_utm, epsg = crs_utm)
land_use_polys_grid    <- lapply(land_use_grid, read_sf_utm, epsg = crs_utm)

# =============================================================================
# 2. COMPUTE COVARIATES (Distance to land use types)
# =============================================================================

cat("\n==========================================================\n")
cat("COMPUTING COVARIATES\n")
cat("==========================================================\n")

# Distance covariates for samples
for (nm in names(land_use_polys_samples)) {
  samples[[paste0("dist_", nm)]] <- min_dist_to(samples, land_use_polys_samples[[nm]])
}

# Distance covariates for prediction grid
for (nm in names(land_use_polys_grid)) {
  pred_grid[[paste0("dist_", nm)]] <- min_dist_to(pred_grid, land_use_polys_grid[[nm]])
}

# Extract coordinates
coords_obs  <- st_coordinates(samples)
coords_pred <- st_coordinates(pred_grid)

# Response variable (log-transformed Pb)
target_col <- intersect(c("Pb_ppm", "Pb__ppm_", "Pb (ppm)"), names(samples))[1]
y_obs <- log1p(st_drop_geometry(samples)[[target_col]])

# Define covariate names
cov_names <- paste0("dist_", names(land_use_polys_samples))

# Standardization function (log1p + scale)
fit_scaler <- function(x) {
  z <- log1p(x)
  list(mu = mean(z, na.rm = TRUE), sd = sd(z, na.rm = TRUE))
}

apply_scaler <- function(x, scaler) {
  (log1p(x) - scaler$mu) / scaler$sd
}

# Fit scalers on training data
scalers <- list()
for (nm in cov_names) {
  scalers[[nm]] <- fit_scaler(st_drop_geometry(samples)[[nm]])
}

# Build design matrices
build_X <- function(sf_data, cov_names, scalers) {
  df <- st_drop_geometry(sf_data)
  X <- matrix(1, nrow = nrow(df), ncol = 1)
  colnames(X) <- "Intercept"
  
  for (nm in cov_names) {
    X <- cbind(X, apply_scaler(df[[nm]], scalers[[nm]]))
  }
  colnames(X)[-1] <- cov_names
  X
}

X_obs  <- build_X(samples, cov_names, scalers)
X_pred <- build_X(pred_grid, cov_names, scalers)

N     <- nrow(X_obs)
P     <- ncol(X_obs)
Npred <- nrow(X_pred)

cat("Observations: N =", N, "\n")
cat("Covariates: P =", P, "\n")
cat("Prediction points: Npred =", Npred, "\n")
cat("Response (log1p Pb): min =", round(min(y_obs), 2), 
    "| max =", round(max(y_obs), 2), 
    "| mean =", round(mean(y_obs), 2), "\n")

# =============================================================================
# 3. EXCEEDANCE THRESHOLD
# =============================================================================

threshold_pb  <- 275  # mg/kg (residential/urban threshold)
threshold_log <- log1p(threshold_pb)

cat("\nExceedance threshold: Pb =", threshold_pb, "mg/kg\n")
cat("Threshold (log scale):", round(threshold_log, 3), "\n")
cat("Observed exceedance rate:", round(mean(y_obs > threshold_log) * 100, 1), "%\n")

# =============================================================================
# 4. MODEL 1: INLA-SPDE
# =============================================================================

cat("\n==========================================================\n")
cat("MODEL 1: INLA-SPDE\n")
cat("==========================================================\n")

t0_inla <- Sys.time()

# Create mesh
bdry_sp <- boundary |> as("Spatial") |> inla.sp2segment()
mesh <- inla.mesh.2d(
  boundary = bdry_sp,
  max.edge = c(1500, 5000),
  offset   = c(1500, 5000),
  cutoff   = 500
)
cat("Mesh nodes:", mesh$n, "\n")

# SPDE model with PC priors
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(5000, 0.5),  # P(range < 5000) = 0.5
  prior.sigma = c(2, 0.5)      # P(sigma > 2) = 0.5
)

# Projection matrices
A_obs  <- inla.spde.make.A(mesh, loc = coords_obs)
A_pred <- inla.spde.make.A(mesh, loc = coords_pred)

s_index <- inla.spde.make.index("s", n.spde = spde$n.spde)

# Data stacks
X_obs_df  <- as.data.frame(X_obs)
X_pred_df <- as.data.frame(X_pred)

stk_obs <- inla.stack(
  data    = list(y = y_obs),
  A       = list(A_obs, 1),
  effects = list(s_index, X_obs_df),
  tag     = "obs"
)

stk_pred <- inla.stack(
  data    = list(y = NA),
  A       = list(A_pred, 1),
  effects = list(s_index, X_pred_df),
  tag     = "pred"
)

stk <- inla.stack(stk_obs, stk_pred)

# Formula
cov_terms <- colnames(X_obs)
formula_inla <- as.formula(
  paste("y ~ 0 +", paste(cov_terms, collapse = " + "), "+ f(s, model = spde)")
)

# Fit model
fit_inla <- inla(
  formula_inla,
  data = inla.stack.data(stk),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
)

t1_inla <- Sys.time()
time_inla <- as.numeric(difftime(t1_inla, t0_inla, units = "secs"))

# Extract predictions
idx_pred <- inla.stack.index(stk, "pred")$data
mu_inla_mean <- fit_inla$summary.fitted.values[idx_pred, "mean"]
mu_inla_sd   <- fit_inla$summary.fitted.values[idx_pred, "sd"]
mu_inla_sd[mu_inla_sd < 1e-6] <- 1e-6

# Exceedance probability (Gaussian approximation)
prob_inla <- 1 - pnorm(threshold_log, mean = mu_inla_mean, sd = mu_inla_sd)

cat("INLA completed in", round(time_inla, 1), "seconds\n")
cat("Exceedance prob: min =", round(min(prob_inla), 4),
    "| max =", round(max(prob_inla), 4),
    "| mean =", round(mean(prob_inla), 4), "\n")

# =============================================================================
# 5. MODEL 2: MCMC (spBayes)
# =============================================================================

cat("\n==========================================================\n")
cat("MODEL 2: MCMC (spBayes)\n")
cat("==========================================================\n")

t0_mcmc <- Sys.time()

# Distance matrix
D_obs <- as.matrix(dist(coords_obs))
maxD  <- max(D_obs)

# Prior specification
phi_lower <- 3 / maxD
phi_upper <- 3 / (0.01 * maxD)

priors_sp <- list(
  beta.Norm   = list(mean = rep(0, P), cov = diag(100, P)),
  sigma.sq.IG = c(2, 1),
  tau.sq.IG   = c(2, 0.5),
  phi.Unif    = c(phi_lower, phi_upper)
)

starting_sp <- list(
  beta     = rep(0, P),
  sigma.sq = var(y_obs) * 0.5,
  tau.sq   = var(y_obs) * 0.5,
  phi      = (phi_lower + phi_upper) / 2
)

tuning_sp <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 0.1)

# Prepare data
dat_sp <- data.frame(y = y_obs, X_obs)
formula_sp <- as.formula(paste("y ~ 0 +", paste(colnames(X_obs), collapse = " + ")))

# Fit model
fit_sp <- spLM(
  formula   = formula_sp,
  coords    = coords_obs,
  data      = dat_sp,
  starting  = starting_sp,
  tuning    = tuning_sp,
  priors    = priors_sp,
  cov.model = "exponential",
  n.samples = 10000,
  verbose   = FALSE
)

# Recover posterior samples
fit_sp_rec <- spRecover(fit_sp, start = 5000, thin = 5, verbose = FALSE)

# Spatial prediction
pred_sp <- spPredict(
  fit_sp_rec,
  pred.covars = X_pred,
  pred.coords = coords_pred,
  verbose     = FALSE
)

t1_mcmc <- Sys.time()
time_mcmc <- as.numeric(difftime(t1_mcmc, t0_mcmc, units = "secs"))

# Extract predictions
y_pred_mcmc  <- pred_sp$p.y.predictive.samples
mu_mcmc_mean <- rowMeans(y_pred_mcmc)
mu_mcmc_sd   <- apply(y_pred_mcmc, 1, sd)

# Exceedance probability (empirical from posterior samples)
prob_mcmc <- rowMeans(y_pred_mcmc > threshold_log)

cat("MCMC completed in", round(time_mcmc, 1), "seconds\n")
cat("Exceedance prob: min =", round(min(prob_mcmc), 4),
    "| max =", round(max(prob_mcmc), 4),
    "| mean =", round(mean(prob_mcmc), 4), "\n")

# =============================================================================
# 6. MODEL 3: VARIATIONAL BAYES (Stan ADVI)
# =============================================================================

cat("\n==========================================================\n")
cat("MODEL 3: VARIATIONAL BAYES (Stan ADVI)\n")
cat("==========================================================\n")

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;
  matrix[N, N] D;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma;
  real<lower=0> rho;
  real<lower=0> tau;
  vector[N] w_raw;
}
transformed parameters {
  vector[N] w;
  vector[N] mu;
  {
    matrix[N, N] K;
    matrix[N, N] L;
    for (i in 1:N) {
      for (j in i:N) {
        K[i, j] = square(sigma) * exp(-D[i, j] / rho);
        K[j, i] = K[i, j];
      }
      K[i, i] += 1e-6;
    }
    L = cholesky_decompose(K);
    w = L * w_raw;
  }
  mu = X * beta + w;
}
model {
  beta ~ normal(0, 10);
  sigma ~ normal(0, 3);
  rho ~ lognormal(log(5000), 1);
  tau ~ normal(0, 2);
  w_raw ~ std_normal();
  y ~ normal(mu, tau);
}
"

t0_vb <- Sys.time()

stan_mod <- stan_model(model_code = stan_code)

fit_vb <- vb(
  stan_mod,
  data = list(N = N, P = P, X = X_obs, y = y_obs, D = D_obs),
  algorithm = "meanfield",
  output_samples = 2000,
  seed = 123
)

t1_vb <- Sys.time()
time_vb <- as.numeric(difftime(t1_vb, t0_vb, units = "secs"))

# Extract posterior means
post <- extract(fit_vb)
beta_vb  <- colMeans(post$beta)
sigma_vb <- mean(post$sigma)
rho_vb   <- mean(post$rho)
tau_vb   <- mean(post$tau)

cat("VB completed in", round(time_vb, 1), "seconds\n")
cat("Estimated params: sigma =", round(sigma_vb, 2), 
    "| rho =", round(rho_vb, 0), 
    "| tau =", round(tau_vb, 2), "\n")

# Kriging prediction with posterior means
D_pn <- rdist(coords_pred, coords_obs)

K_obs <- sigma_vb^2 * exp(-D_obs / rho_vb) + diag(1e-6, N)
K_pn  <- sigma_vb^2 * exp(-D_pn / rho_vb)
Ky    <- K_obs + diag(tau_vb^2, N)

Ky_inv <- solve(Ky)
resid  <- y_obs - X_obs %*% beta_vb
w_pred <- K_pn %*% Ky_inv %*% resid
mu_vb_mean <- as.vector(X_pred %*% beta_vb + w_pred)

# Predictive variance
var_pred <- (sigma_vb^2 + tau_vb^2) - rowSums((K_pn %*% Ky_inv) * K_pn)
var_pred[var_pred < 1e-6] <- 1e-6
mu_vb_sd <- sqrt(var_pred)

# Exceedance probability
prob_vb <- 1 - pnorm(threshold_log, mean = mu_vb_mean, sd = mu_vb_sd)

cat("Exceedance prob: min =", round(min(prob_vb), 4),
    "| max =", round(max(prob_vb), 4),
    "| mean =", round(mean(prob_vb), 4), "\n")

# =============================================================================
# 7. RESULTS SUMMARY
# =============================================================================

cat("\n==========================================================\n")
cat("RESULTS SUMMARY\n")
cat("==========================================================\n")

# Add results to prediction grid
pred_grid$prob_inla <- prob_inla
pred_grid$prob_mcmc <- prob_mcmc
pred_grid$prob_vb   <- prob_vb
pred_grid$diff_inla_mcmc <- prob_inla - prob_mcmc
pred_grid$diff_inla_vb   <- prob_inla - prob_vb

# Summary table
summary_table <- data.frame(
  Method     = c("INLA", "MCMC", "VB"),
  Time_sec   = round(c(time_inla, time_mcmc, time_vb), 1),
  Prob_min   = round(c(min(prob_inla), min(prob_mcmc), min(prob_vb)), 4),
  Prob_max   = round(c(max(prob_inla), max(prob_mcmc), max(prob_vb)), 4),
  Prob_mean  = round(c(mean(prob_inla), mean(prob_mcmc), mean(prob_vb)), 4),
  Pct_gt_10  = round(c(mean(prob_inla > 0.1), mean(prob_mcmc > 0.1), mean(prob_vb > 0.1)) * 100, 1),
  Pct_gt_50  = round(c(mean(prob_inla > 0.5), mean(prob_mcmc > 0.5), mean(prob_vb > 0.5)) * 100, 1)
)

cat("\nComparison Summary:\n")
print(summary_table)

# Correlation between methods
cat("\nCorrelation between methods:\n")
cat("  INLA vs MCMC:", round(cor(prob_inla, prob_mcmc), 4), "\n")
cat("  INLA vs VB:  ", round(cor(prob_inla, prob_vb), 4), "\n")
cat("  MCMC vs VB:  ", round(cor(prob_mcmc, prob_vb), 4), "\n")

# =============================================================================
# 8. VISUALIZATION: HISTOGRAMS
# =============================================================================

cat("\n==========================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("==========================================================\n")

# Histograms of exceedance probabilities
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
hist(prob_inla, breaks = 50, main = "INLA", xlab = "P(Pb > threshold)", 
     col = "steelblue", border = "white")
abline(v = mean(prob_inla), col = "red", lwd = 2, lty = 2)

hist(prob_mcmc, breaks = 50, main = "MCMC", xlab = "P(Pb > threshold)", 
     col = "darkgreen", border = "white")
abline(v = mean(prob_mcmc), col = "red", lwd = 2, lty = 2)

hist(prob_vb, breaks = 50, main = "VB", xlab = "P(Pb > threshold)", 
     col = "darkorange", border = "white")
abline(v = mean(prob_vb), col = "red", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

# =============================================================================
# 9. VISUALIZATION: EXCEEDANCE PROBABILITY MAPS
# =============================================================================

# Map function using geom_sf
make_prob_map <- function(grid_sf, boundary_sf, var, title, subtitle = "") {
  ggplot() +
    geom_sf(data = boundary_sf, fill = "grey90", color = "grey40", linewidth = 0.3) +
    geom_sf(data = grid_sf, aes(color = .data[[var]]), size = 0.6) +
    scale_color_viridis(
      option = "inferno",
      limits = c(0, 1),
      na.value = "grey50",
      name = "P(exceed)"
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text  = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
}

# Generate maps
p_inla <- make_prob_map(
  pred_grid, boundary, "prob_inla",
  paste0("INLA"),
  paste("Threshold:", threshold_pb, "mg/kg")
)

p_mcmc <- make_prob_map(
  pred_grid, boundary, "prob_mcmc",
  paste0("MCMC"),
  paste("Threshold:", threshold_pb, "mg/kg")
)

p_vb <- make_prob_map(
  pred_grid, boundary, "prob_vb",
  paste0("VB"),
  paste("Threshold:", threshold_pb, "mg/kg")
)

# Combine exceedance maps
combined_prob <- plot_grid(p_inla, p_mcmc, p_vb, ncol = 3, align = "h")
print(combined_prob)

ggsave("exceedance_probability_comparison.png", combined_prob, 
       width = 15, height = 5, dpi = 300)
cat("Saved: exceedance_probability_comparison.png\n")

# =============================================================================
# 10. VISUALIZATION: DIFFERENCE MAPS
# =============================================================================

make_diff_map <- function(grid_sf, boundary_sf, var, title) {
  ggplot() +
    geom_sf(data = boundary_sf, fill = "grey90", color = "grey40", linewidth = 0.3) +
    geom_sf(data = grid_sf, aes(color = .data[[var]]), size = 0.6) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-0.5, 0.5),
      oob = squish,
      name = "Δ Prob"
    ) +
    labs(title = title) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text  = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )
}

p_diff1 <- make_diff_map(pred_grid, boundary, "diff_inla_mcmc", "INLA − MCMC")
p_diff2 <- make_diff_map(pred_grid, boundary, "diff_inla_vb", "INLA − VB")

combined_diff <- plot_grid(p_diff1, p_diff2, ncol = 2, align = "h")
print(combined_diff)

ggsave("exceedance_probability_differences.png", combined_diff, 
       width = 10, height = 5, dpi = 300)
cat("Saved: exceedance_probability_differences.png\n")

# =============================================================================
# 11. VISUALIZATION: OBSERVED DATA MAP
# =============================================================================

# Add observed Pb to samples for plotting
samples$Pb_log <- y_obs
samples$exceeds <- factor(ifelse(y_obs > threshold_log, "Yes", "No"), 
                          levels = c("No", "Yes"))

p_obs <- ggplot() +
  geom_sf(data = boundary, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = samples, aes(color = Pb_log, shape = exceeds), size = 2.5) +
  scale_color_viridis(option = "magma", name = "log(Pb+1)") +
  scale_shape_manual(values = c("No" = 16, "Yes" = 17), name = "Exceeds") +
  labs(
    title = "Observed Pb Concentrations",
    subtitle = paste("N =", N, "| Exceedance threshold:", threshold_pb, "mg/kg")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_obs)
ggsave("observed_pb_concentrations.png", p_obs, width = 8, height = 6, dpi = 300)
cat("Saved: observed_pb_concentrations.png\n")

# =============================================================================
# 12. FINAL OUTPUT
# =============================================================================

cat("\n==========================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==========================================================\n")

cat("\nOutput files:\n")
cat("  - exceedance_probability_comparison.png\n")
cat("  - exceedance_probability_differences.png\n")
cat("  - observed_pb_concentrations.png\n")

cat("\nPrediction grid object 'pred_grid' contains:\n")
cat("  - prob_inla: INLA exceedance probabilities\n")
cat("  - prob_mcmc: MCMC exceedance probabilities\n")
cat("  - prob_vb: VB exceedance probabilities\n")
cat("  - diff_inla_mcmc: Difference INLA - MCMC\n")
cat("  - diff_inla_vb: Difference INLA - VB\n")

cat("\nDone!\n")
