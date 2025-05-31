suppressPackageStartupMessages(library(bkmr))
suppressPackageStartupMessages(library(bws))
suppressPackageStartupMessages(library(gWQS))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(qgcomp))
suppressPackageStartupMessages(library(splines))
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
params <- tibble(
  scenario = 1:21,
  n = 500,
  p = 5,
  beta_X = rep(list(
    c(0, 0, 0, 0, 0),          # scenario 1
    c(0.25, 0, 0, 0, 0),          # scenario 2
    c(0.25/4, 0.25/4, 0.25/4, 0.25/4, 0), # scenario 3
    c(0.25, -0.15, 0, 0, 0),     # scenario 4
    c(0, 0, 0, 0, 0),        # scenario 5
    c(0.25, 0, 0, 0, 0),          # scenario 6
    c(0.25, 0.15, 0, 0, 0)         # scenario 7
  ), each = 3),
  beta_C = rep(c(0, 0, 0, 0, 0.25, 0, 0), each = 3),     # scenario 5 has confounding
  beta_X1X1 = rep(c(0, 0, 0, 0, 0, -0.15, 0), each = 3),  # scenario 6 has a nonlinear effect
  beta_X1X2 = rep(c(0, 0, 0, 0, 0, 0, -0.15), each = 3),  # scenario 7 has an interaction effect
  rho_X = rep(c(0, 0.4, 0.7), times = 7), 
  rho_C = rep(c(0, 0, 0, 0, 0.7, 0, 0), each = 3)     # scenario 5 has confounding
)

params <- params[c(7:9, 13:21),]

###############################################################
## start simulation code
###############################################################
nsim <- 100
# define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 21
  nsim = 1
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

# define simulation scenario
#for (scenario in 11){

param <- params[scenario, ]

# generate a random seed for each simulated dataset
seed <- floor(runif(nsim, 1, 900))
results = vector("list", length = nsim)

# run simulations
for(i in 1:nsim){
  cat("scenario:", scenario, ", i:", i, "\n")
  set.seed(seed[i])
  
  ####################
  # simulate data
  simdata <- simulate_data(n = param$n,
                           p = param$p,
                         beta_X = param$beta_X[[1]],
                         beta_C = param$beta_C,
                         beta_X1X1 = param$beta_X1X1,
                         beta_X1X2 = param$beta_X1X2,
                         rho_X = param$rho_X,
                         rho_C = param$rho_C)
  
  ####################
  # WQS with a positive indice
  fit.wqs <- tryCatch(
    {
      gwqs(y ~ wqs, mix_name = paste0("X", 1:5), data = simdata,
           q = 4, validation = 0.6, b1_pos = TRUE, b = 100, rh = 100,
           family = "gaussian")
    },
    error = function(e) {
      message(paste("Simulation", i, "failed with error:", e$message))
      return(NULL)
    }
  )
  
  if (is.null(fit.wqs)) next
  
  # WQS with two indices
  fit.wqs2 <- tryCatch(
    {
      gwqs(y ~ pwqs + nwqs, mix_name = paste0("X", 1:5), data = simdata, 
           q = 4, validation = 0.6, b1_pos = TRUE, b = 100, rh = 100,
           family = "gaussian")
    },
    error = function(e) {
      message(paste("Simulation", i, "failed with error:", e$message))
      return(NULL)
    }
  )
  
  if (is.null(fit.wqs2)) next
  
  # qgcomp.noboot
  fit.qgcomp <- qgcomp.glm.noboot(y ~ ., dat = simdata, 
                                  family = gaussian(), q = 4, bayes = TRUE)
  
  
  # qgcomp.boot
  if (param$beta_X1X1 != 0){
    fit.qgcomp.boot <- qgcomp.glm.boot(y ~ bs(X1) + X2 + X3 + X4 + X5, dat = simdata, 
                                       expnms = paste0("X", 1:5),
                                       family = gaussian(), q = 4, B = 200, degree = 2)
  } else if (param$beta_X1X2 != 0){
    fit.qgcomp.boot <- qgcomp.glm.boot(y ~ bs(X1) * bs(X2) + X3 + X4 + X5, dat = simdata, 
                                       expnms = paste0("X", 1:5),
                                       family = gaussian(), q = 4, B = 200, degree = 2)
  }
  
  # BKMR
  fit.bkmr <- kmbayes(y = simdata$y, Z = simdata[, -1], 
                      family = "gaussian", iter = 2000, verbose = FALSE, varsel = TRUE)
  
  # BWS
  fit.bws <- bws(iter = 2000, y = simdata$y, X = simdata[, -1], family = "gaussian")
  
  ####################
  # extract mixture effect and individual weights estimates
  res.wqs <- extract_estimates(model = fit.wqs, method = "WQS")
  res.wqs2 <- extract_estimates(model = fit.wqs2, method = "WQS2")
  res.qgcomp <- extract_estimates(model = fit.qgcomp, method = "qgcomp.noboot")
  if (param$beta_X1X1 != 0 | param$beta_X1X2 != 0){
    res.qgcomp.boot <- extract_estimates(model = fit.qgcomp.boot, method = "qgcomp.boot")
  }
  res.bkmr <- extract_estimates(model = fit.bkmr, method = "BKMR")
  res.bws <- extract_estimates(model = fit.bws, method = "BWS")
  
  if (param$beta_X1X1 != 0 | param$beta_X1X2 != 0){
    df_coef <- bind_rows(res.wqs[[1]], res.wqs2[[1]], res.qgcomp[[1]], res.qgcomp.boot[[1]], res.bws[[1]]) %>%
      bind_cols(param, seed = seed[i])
    df_coef$AIC_qgcomp = AIC(fit.qgcomp)
    df_coef$AIC_qgcomp_boot = AIC(fit.qgcomp.boot)
  } else {
    df_coef <- bind_rows(res.wqs[[1]], res.wqs2[[1]], res.qgcomp[[1]], res.bws[[1]]) %>%
      bind_cols(param, seed = seed[i])
  }
  #rownames(df_coef) <- NULL
  df_weights <- bind_rows(res.wqs[[2]], res.wqs2[[2]], res.qgcomp[[2]], res.bkmr[[2]], res.bws[[2]]) %>%
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
date = gsub("-", "", Sys.Date())
dir.create(file.path(here("results"), date), showWarnings = FALSE)

filename = paste0(here("results", date), "/", scenario, ".RDA")
save(results, file = filename)
#}
