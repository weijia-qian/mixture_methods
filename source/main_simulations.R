suppressPackageStartupMessages(library(bkmr))
suppressPackageStartupMessages(library(bws))
suppressPackageStartupMessages(library(gWQS))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(qgcomp))
suppressPackageStartupMessages(library(tidyverse))

wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################

source(here("source", "simulate_data.R"))
source(here("source", "extract_estimates.R"))

###############################################################
## set simulation design elements
###############################################################

nsim = 10

params <- tibble(
  n = 500,
  beta_X = list(
    c(0, 0, 0, 0, 0),          # scenario 1
    c(1, 0, 0, 0, 0),          # scenario 2
    c(0.25, 0.25, 0.25, 0.25, 0), # scenario 3
    c(0.5, 0.5, 0, 0, 0),      # scenario 4
    c(0.5, -0.5, 0, 0, 0),     # scenario 5
    c(0.5, 0, 0, 0, 0),        # scenario 6
    c(1, 0, 0, 0, 0),          # scenario 7
    c(1, 0.5, 0, 0, 0)         # scenario 8
  ),
  beta_C = c(0, 0, 0, 0, 0, 0.5, 0, 0),     # scenario 6 has confounding
  beta_X1X1 = c(0, 0, 0, 0, 0, 0, 0.2, 0),  # scenario 7 has a nonlinear effect
  beta_X1X2 = c(0, 0, 0, 0, 0, 0, 0, 0.2),  # scenario 8 has an interaction effect
  rho_X = c(0, 0, 0, 0.25, 0.25, 0, 0, 0),  # scenarios 4 and 5 have correlated exposures
  rho_C = c(0, 0, 0, 0, 0, 0.25, 0, 0)      # scenario 6 has a correlated confounder
)

###############################################################
## start simulation code
###############################################################

# define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 1
  nsim = 2
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

# define simulation scenario
param <- params[scenario, ]

# generate a random seed for each simulated dataset
seed <- floor(runif(nsim, 1, 900))
results = vector("list", length = nsim)

# run simulations
for(i in 1:nsim){
  
  set.seed(seed[i])
  
  ####################
  # simulate data
  simdata <- simulate_data(n = param$n,
                         beta_X = param$beta_X[[1]],
                         beta_C = param$beta_C,
                         beta_X1X1 = param$beta_X1X1,
                         beta_X1X2 = param$beta_X1X2,
                         rho_X = param$rho_X,
                         rho_C = param$rho_C)
  
  ####################
  # WQS
  fit.wqs <- gwqs(y ~ wqs, mix_name = paste0("X", 1:5), data = simdata, 
                  q = NULL, validation = 0.6, b1_pos = TRUE, b = 100, rh = 100,
                  family = "gaussian")
  
  # qgcomp.noboot
  fit.qgcomp <- qgcomp.glm.noboot(y ~ X1 + X2 + X3 + X4 + X5, dat = simdata, 
                                  family = gaussian(), q = NULL, bayes = TRUE)
  
  # qgcomp.boot
  fit.qgcomp.boot <- qgcomp.glm.boot(y ~ X1 + X2 + X3 + X4 + X5, dat = simdata, 
                                     family = gaussian(), q = NULL, bayes = TRUE, B = 200)
  
  # BKMR
  fit.bkmr <- kmbayes(y = simdata$y, Z = simdata[, -1], 
                      family = "gaussian", iter = 1000, verbose = FALSE, varsel = TRUE)
  
  # BWS
  fit.bws <- bws(iter = 1000, y = simdata$y, X = simdata[, -1], family = "gaussian")
  
  ####################
  # extract mixture effect and individual weights estimates
  res.wqs <- extract_estimates(model = fit.wqs, method = "WQS")
  res.qgcomp <- extract_estimates(model = fit.qgcomp, method = "qgcomp.noboot")
  res.qgcomp.boot <- extract_estimates(model = fit.qgcomp.boot, method = "qgcomp.boot")
  res.bkmr <- extract_estimates(model = fit.bkmr, method = "BKMR")
  res.bws <- extract_estimates(model = fit.bws, method = "BWS")
  
  df_coef <- bind_rows(res.wqs[[1]], res.qgcomp[[1]], res.qgcomp.boot[[1]], res.bws[[1]]) %>%
    bind_cols(param, seed = seed[i])
  rownames(df_coef) <- NULL
  df_weights <- bind_rows(res.wqs[[2]], res.qgcomp[[2]], res.bkmr[[2]], res.bws[[2]]) %>%
    bind_cols(param, seed = seed[i])
  rownames(df_weights) <- NULL
  df_bkmr <- res.bkmr[[1]]  %>%
    bind_cols(param, seed = seed[i])
  rownames(df_bkmr) <- NULL
  
  ####################
  # store results
  res <- list(coef = df_coef, weights = df_weights, bkmr = df_bkmr)
  
  results[[i]] = res

} # end for loop

####################
# save results
filename <- paste0("scenario_", scenario, ".RDA")
save(results, file = here("results", filename))


