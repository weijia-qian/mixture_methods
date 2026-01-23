extract_estimates <- function(model, method, simdata = NULL) {
  
  if (method == "WQS") {
    
    # mixture effect
    coef_all <- coef(model)
    term <- names(coef_all)[2]
    est <- unname(coef_all[2])
    
    ci <- suppressWarnings(confint(model))
    lb <- unname(ci[2, 1])
    ub <- unname(ci[2, 2])
    
    tmp_coef <- data.frame(
      method = method,
      term = term,
      estimate = est,
      lb = lb,
      ub = ub
    )
    
    weights <- model$final_weights
    tmp_weights <- data.frame(
      method = method,
      mix_name = weights$mix_name,
      estimate = weights$Estimate,
      lb = weights$`2.5%`,
      ub = weights$`97.5%`
    )
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "WQS2") {
    
    coef_all <- coef(model)
    terms <- names(coef_all)[-1]
    est <- unname(coef_all[-1])
    
    ci <- suppressWarnings(confint(model))
    lb <- unname(ci[-1, 1])
    ub <- unname(ci[-1, 2])
    
    tmp_coef <- data.frame(
      method = c("WQS2: pos", "WQS2: neg"),
      term = terms,
      estimate = est,
      lb = lb,
      ub = ub
    )
    
    tmp_weights <- rbind(
      data.frame(
        method = "WQS2: pos",
        mix_name = model$final_weights$mix_name,
        estimate = model$final_weights$`Estimate pos`,
        lb = model$final_weights$`2.5% pos`,
        ub = model$final_weights$`97.5% pos`
      ),
      data.frame(
        method = "WQS2: neg",
        mix_name = model$final_weights$mix_name,
        estimate = model$final_weights$`Estimate neg`,
        lb = model$final_weights$`2.5% neg`,
        ub = model$final_weights$`97.5% neg`
      )
    )
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "qgcomp.noboot") {
    
    # mixture effect (psi1)
    coef_all <- coef(model)
    term <- names(coef_all)[2]
    est <- unname(coef_all[2])
    
    ci <- suppressWarnings(confint(model))
    lb <- unname(ci[2, 1])
    ub <- unname(ci[2, 2])
    
    tmp_coef <- data.frame(
      method = method,
      term = term,
      estimate = est,
      lb = lb,
      ub = ub
    )
    
    # "weights" (signed contribution proxy)
    # NOTE: qgcomp weights are typically interpreted as directional contributions.
    pos_weights_names <- names(model$pos.weights)
    pos_weights <- model$pos.weights * model$pos.psi / (model$pos.size + model$neg.size)
    
    neg_weights_names <- names(model$neg.weights)
    neg_weights <- model$neg.weights * model$neg.psi / (model$pos.size + model$neg.size)
    
    tmp_weights <- data.frame(
      method = method,
      mix_name = c(pos_weights_names, neg_weights_names),
      estimate = c(pos_weights, neg_weights),
      lb = NA_real_,
      ub = NA_real_
    )
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "qgcomp.boot") {
    
    # Boot object: store named parameters so downstream code can pick psi1, psi2, etc.
    coef_all <- coef(model)
    ci <- suppressWarnings(confint(model))
    
    tmp_coef <- data.frame(
      method = method,
      term = names(coef_all),
      estimate = unname(coef_all),
      lb = unname(ci[, 1]),
      ub = unname(ci[, 2])
    )
    
    tmp_weights <- data.frame(
      method = method,
      mix_name = NA_character_,
      estimate = NA_real_,
      lb = NA_real_,
      ub = NA_real_
    )
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "BKMR") {
    
    stopifnot(!is.null(simdata))
    
    weights <- ExtractPIPs(model)
    tmp_weights <- data.frame(
      method = method,
      mix_name = weights$variable,
      estimate = weights$PIP,
      lb = NA_real_,
      ub = NA_real_
    )
    
    risks.overall <- OverallRiskSummaries(
      fit = model,
      y = simdata$y,
      Z = simdata[, -1],
      qs = seq(0, 1, by = 0.05),
      q.fixed = 0.5,
      method = "exact"
    )
    
    tmp_bkmr <- data.frame(
      quantile = risks.overall$quantile,
      est = risks.overall$est,
      sd = risks.overall$sd
    )
    
    return(list(tmp_bkmr, tmp_weights))
    
  } else if (method == "BWS") {
    
    df_bws <- as.data.frame(model, pars = c("theta1", "w"))
    
    # mixture effect
    est <- mean(df_bws$theta1)
    lb <- unname(quantile(df_bws$theta1, 0.025))
    ub <- unname(quantile(df_bws$theta1, 0.975))
    
    tmp_coef <- data.frame(
      method = method,
      term = "theta1",
      estimate = est,
      lb = lb,
      ub = ub
    )
    
    weights <- colMeans(df_bws)[-1]
    weights_lb <- apply(df_bws[, -1, drop = FALSE], 2, quantile, probs = 0.025)
    weights_ub <- apply(df_bws[, -1, drop = FALSE], 2, quantile, probs = 0.975)
    p <- length(weights)
    
    tmp_weights <- data.frame(
      method = method,
      mix_name = paste0("X", 1:p),
      estimate = unname(weights),
      lb = unname(weights_lb),
      ub = unname(weights_ub)
    )
    
    return(list(tmp_coef, tmp_weights))
  }
  
  stop("Unknown method: ", method)
}
