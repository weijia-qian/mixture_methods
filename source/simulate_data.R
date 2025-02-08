###############################################################
## function to simulate data
###############################################################
simulate_data <- function(
    n = 500,                      # sample size
    beta_0 = 0,                       # model intercept
    beta_X = c(0, 0, 0, 0, 0),    # beta coefficients for X
    beta_X1X1 = 0,              # beta coefficient for the X1X2 interaction
    beta_X1X2 = 0,              # beta coefficient for the X1X2 interaction
    beta_C = 0,                   # beta coefficient for the unmeasured confounder (C)
    rho_X1X2 = 0,                    # correlation between X1 and X2
    rho_C = 0.75                 # correlation between X1 and C
){
  
  p = length(beta_X) # number of predictors

    # Simulate X
    X = matrix(nrow = n, ncol = p) 
    for (k in 1:p) {
      X[, k] <- rnorm(n)
    }
    colnames(X) = paste0("X", 1:p)
    
    # Introduce correlation between X1 and X2
    if (rho_X1X2 != 0){
      X[, 2] <- (rho_X1X2 * X[, 1]) + sqrt(1 - rho_X1X2 * rho_X1X2) * rnorm(n)
    }
    
    # Simulate the unmeasured confounder C
    C <- (rho_C * X[, 1]) + sqrt(1 - rho_C * rho_C) * rnorm(n)
    
    # Simulate outcome
    y <- beta_0 + X %*% beta_X + X[, 1]^2 * beta_X1X1 + (X[, 1] * X[, 2]) * beta_X1X2 + C * beta_C + rnorm(n)
    
    res <- data.frame(y = y, 
                      X)
  return(res)
}
