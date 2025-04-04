###############################################################
## function to simulate data
###############################################################
simulate_data <- function(
    n = 500,                      # sample size
    beta_0 = 0,                   # model intercept
    beta_X = c(0, 0, 0, 0, 0),    # beta coefficients for X
    beta_X1X1 = 0,                # beta coefficient for the X1X2 interaction
    beta_X1X2 = 0,                # beta coefficient for the X1X2 interaction
    beta_C = 0,                   # beta coefficient for the unmeasured confounder (C)
    rho_X = 0,                    # pairwise correlation between X's
    rho_C = 0.7                # correlation between X1 and C
){
    p = length(beta_X) # number of predictors

    # Simulate X
    Sigma <- matrix(rho_X, nrow = p, ncol = p)
    diag(Sigma) <- 1
    
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    colnames(X) = paste0("X", 1:p)
    
    # Simulate the unmeasured confounder C
    C <- (rho_C * X[, 1]) + sqrt(1 - rho_C * rho_C) * rnorm(n)
    
    # Simulate outcome
    y <- beta_0 + X %*% beta_X + X[, 1]^2 * beta_X1X1 + (X[, 1] * X[, 2]) * beta_X1X2 + C * beta_C + rnorm(n)
    
    res <- data.frame(y = y, 
                      X)
  return(res)
}