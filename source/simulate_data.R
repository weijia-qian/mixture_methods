###############################################################
## function to simulate data
###############################################################
simulate_data <- function(
    n = 500,
    p = 5,
    scenario = c("null", "single", "homogeneous", "heterogeneous",
                 "nonlinear", "interactive"),
    rho_X = 0,
    sigma = 1,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  scenario <- match.arg(scenario)
  
  ## ---------------------------
  ## 1. Simulate exposures
  ## ---------------------------
  Sigma <- matrix(rho_X, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  
  q <- qnorm(0.75)  # 75th percentile
  
  ## ---------------------------
  ## 2. Generate outcome
  ## ---------------------------
  y <- numeric(n)
  
  if (scenario == "null") {
    
    y <- rnorm(n, 0, sigma)
    
  }
  
  if (scenario == "single") {
    
    # One active exposure, effect size = 0.25
    beta <- c(0.25, rep(0, p - 1))
    y <- X %*% beta + rnorm(n, 0, sigma)
    
  }
  
  if (scenario == "homogeneous") {
    
    # Four equal positive effects summing to 0.25
    beta <- c(rep(0.0625, 4), rep(0, p - 4))
    y <- X %*% beta + rnorm(n, 0, sigma)
    
  }
  
  if (scenario == "heterogeneous") {
    
    # Mixed directions, net effect = 0.10
    beta <- c(0.25, -0.15, rep(0, p - 2))
    y <- X %*% beta + rnorm(n, 0, sigma)
    
  }
  
  if (scenario == "nonlinear") {
    
    # Smooth additive nonlinearity (not tailored to any method)
    y <- 0.35 * X[, 1] -
      0.15 * X[, 1]^2 +
      0.20 * sin(pi * X[, 2]) +
      rnorm(n, 0, sigma)
  }
  
  if (scenario == "interactive") {
    
    # Strong threshold-based interaction
    I1 <- as.numeric(X[, 1] > q)
    I2 <- as.numeric(X[, 2] > q)
    
    y <- 0.25 * X[, 1] * I1 +
      0.15 * X[, 2] * I2 -
      0.15 * X[, 1] * X[, 2] * I1 * I2 +
      rnorm(n, 0, sigma)
  }
  
  ## ---------------------------
  ## 3. Return data
  ## ---------------------------
  data.frame(y = y, X)
}

