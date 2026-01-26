suppressPackageStartupMessages(library(bkmr))
suppressPackageStartupMessages(library(bws))
suppressPackageStartupMessages(library(gWQS))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(qgcomp))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(tidyverse))

# Optional: avoid future warnings from gWQS nested parallelism
suppressPackageStartupMessages(library(future))
future::plan(future::sequential)

wd <- getwd()
doLocal <- substring(wd, 2, 6) == "Users"

###############################################################
## define or source functions used in code below
###############################################################
source(here("source", "simulate_data.R"))
source(here("source", "extract_estimates.R"))

###############################################################
## set simulation design elements
###############################################################
scenarios <- c("single", "homogeneous", "heterogeneous")
rho_levels <- c(0, 0.4, 0.7)
sigma_levels <- c(1.0)
p_levels <- c(10, 20)

params <- tidyr::crossing(
  scenario = scenarios,
  rho_X = rho_levels,
  sigma = sigma_levels,
  n = 500,
  p = p_levels
) %>%
  dplyr::mutate(batch = dplyr::row_number()) %>%
  dplyr::select(batch, dplyr::everything())

pick_top_exposures <- function(dat, Xnms, K = 4) {
  # simple, common screening: absolute correlation with outcome
  # (fast + works even when qgcomp weights are unstable)
  cors <- sapply(Xnms, function(x) suppressWarnings(cor(dat[[x]], dat$y)))
  cors[is.na(cors)] <- 0
  names(sort(abs(cors), decreasing = TRUE))[seq_len(min(K, length(Xnms)))]
}

###############################################################
## start simulation code
###############################################################
nsim <- 25

if (doLocal) {
  batch <- 1
  nsim <- 1
} else {
  batch <- as.numeric(commandArgs(trailingOnly = TRUE))
}

param <- params[batch, ]

# reproducible seeds per batch
set.seed(10000 + batch)
seed <- sample.int(1e8, nsim)

results <- vector("list", length = nsim)
time <- vector("numeric", length = nsim)
Xnms <- paste0("X", 1:param$p)

for (i in 1:nsim) {
  cat("batch:", batch, ", i:", i, "\n")
  tic()
  ####################
  # simulate data
  simdata <- simulate_data(
    n = param$n,
    p = param$p,
    scenario = param$scenario,
    rho_X = param$rho_X,
    sigma = param$sigma,
    seed = seed[i]
  )
  
  ####################
  # WQS with a positive index
  fit.wqs <- tryCatch(
    {
      gwqs(
        y ~ wqs,
        mix_name = Xnms,
        data = simdata,
        q = 4,
        validation = 0.6,
        b1_pos = TRUE,
        b = 100,
        rh = 100,
        family = "gaussian"
      )
    },
    error = function(e) {
      message(paste("WQS failed (batch", batch, "sim", i, "):", e$message))
      return(NULL)
    }
  )
  if (is.null(fit.wqs)) next
  
  ####################
  # WQS with two indices
  fit.wqs2 <- tryCatch(
    {
      gwqs(
        y ~ pwqs + nwqs,
        mix_name = Xnms,
        data = simdata,
        q = 4,
        validation = 0.6,
        b1_pos = TRUE,
        b = 100,
        rh = 100,
        family = "gaussian"
      )
    },
    error = function(e) {
      message(paste("WQS2 failed (batch", batch, "sim", i, "):", e$message))
      return(NULL)
    }
  )
  if (is.null(fit.wqs2)) next
  
  ####################
  # qgcomp baseline (linear)
  fit.qgcomp <- tryCatch(
    {
      qgcomp.glm.noboot(
        y ~ .,
        dat = simdata,
        family = gaussian(),
        q = 4,
        bayes = TRUE
      )
    },
    error = function(e) {
      message(paste("qgcomp.noboot failed (batch", batch, "sim", i, "):", e$message))
      return(NULL)
    }
  )
  if (is.null(fit.qgcomp)) next
  
  ####################
  # qgcomp extended (practice-style)
  # - quadratic terms for ALL exposures
  # - for interactive scenario: interactions among TOP K exposures (screened)
  fit.qgcomp.ext <- NULL
  if (param$scenario %in% c("nonlinear", "interactive")) {
    
    Xnms <- paste0("X", 1:param$p)
    quad_terms <- paste0("I(", Xnms, "^2)", collapse = " + ")
    
    if (param$scenario == "nonlinear") {
      
      f_ext <- stats::as.formula(
        paste0("y ~ ", paste(Xnms, collapse = " + "), " + ", quad_terms)
      )
      
    } else {
      # interactive: choose top K exposures then include only their pairwise interactions
      topK <- pick_top_exposures(simdata, Xnms, K = 4)
      int_terms <- paste0("(", paste(topK, collapse = " + "), ")^2")
      
      # include full linear terms for all X + quadratics for all X + limited interactions among topK
      f_ext <- stats::as.formula(
        paste0("y ~ ", paste(Xnms, collapse = " + "),
               " + ", quad_terms,
               " + ", int_terms)
      )
    }
    
    fit.qgcomp.ext <- tryCatch(
      {
        qgcomp.glm.boot(
          f_ext,
          dat = simdata,
          expnms = Xnms,
          family = gaussian(),
          q = 4,
          B = 100
        )
      },
      error = function(e) {
        message(paste("qgcomp.boot failed (batch", batch, "sim", i, "):", e$message))
        return(NULL)
      }
    )
  }
  
  ####################
  # BKMR
  fit.bkmr <- tryCatch(
    {
      kmbayes(
        y = simdata$y,
        Z = simdata[, -1],
        family = "gaussian",
        iter = 2000,
        verbose = FALSE,
        varsel = TRUE
      )
    },
    error = function(e) {
      message(paste("BKMR failed (batch", batch, "sim", i, "):", e$message))
      return(NULL)
    }
  )
  if (is.null(fit.bkmr)) next
  
  ####################
  # BWS
  fit.bws <- tryCatch(
    {
      bws::bws(
        iter = 2000,                 # use 2000 to keep runtime sane; change if you want
        y = simdata$y,
        X = simdata[, -1],
        family = "gaussian",
        refresh = 0                 # suppress Stan progress printing
      )
    },
    error = function(e) {
      message(paste("BWS failed (batch", batch, "sim", i, "):", e$message))
      return(NULL)
    }
  )
  if (is.null(fit.bws)) next
  
  ####################
  # extract estimates
  res.wqs   <- extract_estimates(model = fit.wqs, method = "WQS")
  res.wqs2  <- extract_estimates(model = fit.wqs2, method = "WQS2")
  res.qg    <- extract_estimates(model = fit.qgcomp, method = "qgcomp.noboot")
  res.bkmr  <- extract_estimates(model = fit.bkmr, method = "BKMR", simdata = simdata)
  res.bws   <- extract_estimates(model = fit.bws, method = "BWS")
  
  res.qgext <- NULL
  if (!is.null(fit.qgcomp.ext)) {
    res.qgext <- extract_estimates(model = fit.qgcomp.ext, method = "qgcomp.boot")
  }
  
  ####################
  # coefficients
  df_coef <- dplyr::bind_rows(
    res.wqs[[1]],
    res.wqs2[[1]],
    res.qg[[1]],
    if (!is.null(res.qgext)) res.qgext[[1]],
    res.bws[[1]]
  ) %>%
    dplyr::mutate(
      batch = param$batch,
      scenario = param$scenario,
      n = param$n,
      p = param$p,
      rho_X = param$rho_X,
      sigma = param$sigma,
      seed = seed[i]
    )
  
  # attach BWS diagnostics to all rows (handy for filtering later)
  bws_diag <- attr(fit.bws, "bws_diag")
  if (!is.null(bws_diag)) {
    df_coef$bws_rhat_max <- bws_diag$rhat_max
    df_coef$bws_ess_min  <- bws_diag$ess_min
    df_coef$bws_refit    <- bws_diag$refit
    df_coef$bws_ok       <- bws_diag$ok
  }
  
  # store AIC comparison if extended qgcomp exists
  if (!is.null(fit.qgcomp.ext)) {
    df_coef$AIC_qgcomp     <- suppressWarnings(AIC(fit.qgcomp))
    df_coef$AIC_qgcomp_ext <- suppressWarnings(AIC(fit.qgcomp.ext))
  }
  
  ####################
  # weights/PIPs
  df_weights <- dplyr::bind_rows(
    res.wqs[[2]],
    res.wqs2[[2]],
    res.qg[[2]],
    res.bkmr[[2]],
    res.bws[[2]]
  ) %>%
    dplyr::mutate(
      batch = param$batch,
      scenario = param$scenario,
      n = param$n,
      p = param$p,
      rho_X = param$rho_X,
      sigma = param$sigma,
      seed = seed[i]
    )
  rownames(df_weights) <- NULL
  
  ####################
  # BKMR overall risk summaries
  df_bkmr <- res.bkmr[[1]] %>%
    dplyr::mutate(
      batch = param$batch,
      scenario = param$scenario,
      n = param$n,
      p = param$p,
      rho_X = param$rho_X,
      sigma = param$sigma,
      seed = seed[i]
    )
  rownames(df_bkmr) <- NULL
  
  results[[i]] <- list(coef = df_coef, weights = df_weights, bkmr = df_bkmr)
  time_stamp <- toc(quiet = TRUE)
  time[i] <- time_stamp$toc - time_stamp$tic
}

####################
# save results
# date <- gsub("-", "", Sys.Date())
dir.create(file.path(here("results"), "k10_20"), showWarnings = FALSE, recursive = TRUE)

filename <- file.path(here("results", "k10_20"), paste0(batch, "_k10_20_batch2.RDA"))
save(results, file = filename)
