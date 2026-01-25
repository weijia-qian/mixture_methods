############################################################
# Monte Carlo signal variance, SNR, and R^2 table (Table S1)
############################################################

suppressPackageStartupMessages({
  library(MASS)      # mvrnorm
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# --------------------------
# Settings
# --------------------------
set.seed(42)
N_mc   <- 200000
p_exp  <- 5
rhos   <- c(0, 0.4, 0.7)
sigmas <- c(0.5, 1, 2)

# interactive threshold (must match simulate_data.R)
q75 <- qnorm(0.75)

# --------------------------
# Draw correlated exposures
# --------------------------
sim_X <- function(n, p, rho) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  X
}

# --------------------------
# Signal function f(X) by scenario
# (Update these to match your simulate_data.R exactly)
# --------------------------
signal_f <- function(X, scenario) {
  if (scenario == "single") {
    return(0.25 * X[, 1])
  }
  if (scenario == "homogeneous") {
    return(0.0625 * rowSums(X[, 1:4, drop = FALSE]))
  }
  if (scenario == "heterogeneous") {
    return(0.25 * X[, 1] - 0.15 * X[, 2])
  }
  if (scenario == "nonlinear") {
    return(0.35 * X[, 1] - 0.15 * (X[, 1]^2) + 0.20 * sin(pi * X[, 2]))
  }
  if (scenario == "interactive") {
    I1 <- as.numeric(X[, 1] >= q75)
    I2 <- as.numeric(X[, 2] >= q75)
    return(0.25 * I1 + 0.15 * I2 - 0.15 * I1 * I2)
  }
  stop("Unknown scenario: ", scenario)
}

# --------------------------
# Compute Var(signal), SNR, R^2 for each (scenario, rho)
# --------------------------
calc_one <- function(scenario, rho, n_mc, p_dim, sigmas_vec) {
  X  <- sim_X(n = n_mc, p = p_dim, rho = rho)
  fX <- signal_f(X, scenario)
  
  v_signal <- var(fX)
  
  snr_vals <- v_signal / (sigmas_vec^2)
  r2_vals  <- v_signal / (v_signal + (sigmas_vec^2))
  
  tibble(
    Scenario = str_to_title(scenario),
    rho = rho,
    `Var(f(X))` = v_signal,
    SNR = paste(sprintf("%.3f", snr_vals), collapse = " / "),
    `R^2` = paste(sprintf("%.3f", r2_vals), collapse = " / ")
  )
}

scenarios <- c("single", "homogeneous", "heterogeneous", "nonlinear", "interactive")

# Use tidyr::expand_grid (recommended) + purrr::pmap
res_tbl <- tidyr::expand_grid(
  scenario = scenarios,
  rho = rhos
) %>%
  mutate(
    out = purrr::map2(
      scenario, rho,
      ~ calc_one(.x, .y, n_mc = N_mc, p_dim = p_exp, sigmas_vec = sigmas)
    )
  ) %>%
  select(out) %>%
  tidyr::unnest(out) %>%
  arrange(factor(Scenario, levels = str_to_title(scenarios)), rho) %>%
  mutate(`Var(f(X))` = round(`Var(f(X))`, 3))

print(res_tbl)

# Optional: save
# write.csv(res_tbl, "TableS1_signal_SNR_R2.csv", row.names = FALSE)
