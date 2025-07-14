extract_estimates <- function(model, method){
  if (method == "WQS"){
    # extract estimated mixture effect
    coef <- coef(model)[2]
    coef_lb <- confint(model)[2, 1]
    coef_ub <- confint(model)[2, 2]
    
    # extract estimated weights
    weights <- model$final_weights
    
    # save results
    tmp_coef <- data.frame(method = method,
                           estimate = coef,
                           lb = coef_lb,
                           ub = coef_ub)
      
    tmp_weights <- data.frame(method = method,
                                mix_name = weights$mix_name,
                                estimate = weights$Estimate,
                                lb = weights$`2.5%`,
                                ub = weights$`97.5%`)
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "WQS2"){
    # extract estimated mixture effect
    coef <- coef(model)[-1]
    coef_lb <- confint(model)[-1, 1]
    coef_ub <- confint(model)[-1, 2]
    
    # extract estimated weights
    weights_pos <- data.frame(method = paste0(method, ": pos"),
                                mix_name = model$final_weights$mix_name,
                                estimate = model$final_weights$`Estimate pos`,
                                lb = model$final_weights$`2.5% pos`,
                                ub = model$final_weights$`97.5% pos`)
    weights_neg <- data.frame(method = paste0(method, ": neg"),
                                mix_name = model$final_weights$mix_name,
                                estimate = model$final_weights$`Estimate neg`,
                                lb = model$final_weights$`2.5% neg`,
                                ub = model$final_weights$`97.5% neg`)
      
    # save results
    tmp_coef <- data.frame(method = c("WQS2: pos", "WQS2: neg"),
                           estimate = coef,
                           lb = coef_lb,
                           ub = coef_ub)
    
    tmp_weights <- rbind(weights_pos, weights_neg)
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "qgcomp.noboot"){
    # extract estimated mixture effect
    coef <- coef(model)[2]
    coef_lb <- confint(model)[2, 1]
    coef_ub <- confint(model)[2, 2]
    
    # extract estimated weights (normalized by effect size)
    pos_weights_names <- names(model$pos.weights)
    pos_weights <- model$pos.weights * model$pos.psi / (model$pos.size + model$neg.size)
    neg_weights_names <- names(model$neg.weights)
    neg_weights <- model$neg.weights * model$neg.psi / (model$pos.size + model$neg.size)
    
    # save results
    tmp_coef <- data.frame(method = method,
                           estimate = coef,
                           lb = coef_lb,
                           ub = coef_ub)
    
    tmp_weights <- data.frame(method = method,
                              mix_name = c(pos_weights_names, neg_weights_names),
                              estimate = c(pos_weights, neg_weights),
                              lb = NA, ub = NA)
    
    return(list(tmp_coef, tmp_weights))
    
  } else if (method == "qgcomp.boot"){
    # extract estimated mixture effect
    coef <- coef(model)[-1]
    coef_lb <- confint(model)[-1, 1]
    coef_ub <- confint(model)[-1, 2]
    
    # save results
    tmp_coef <- data.frame(method = method,
                           estimate = coef,
                           lb = coef_lb,
                           ub = coef_ub)
    
    return(list(tmp_coef))
    
  } else if (method == "BKMR"){
    # estimated posterior inclusion probabilities
    weights <- ExtractPIPs(model)
    tmp_weights <- data.frame(method = method,
                              mix_name = weights$variable,
                              estimate = weights$PIP,
                              lb = NA, ub = NA)
  
    # summary statistics of the predictor-response function
    risks.overall <- OverallRiskSummaries(fit = model, y = simdata$y, Z = simdata[, -1], 
                                          qs = seq(0, 1, by = 0.05), 
                                          q.fixed = 0.5, method = "exact")
    tmp_bkmr <- data.frame(quantile = risks.overall$quantile,
                           est = risks.overall$est,
                           sd = risks.overall$sd)
    
    return(list(tmp_bkmr, tmp_weights))
    
  } else if (method == "BWS"){
    # extract estimated mixture effect
    df_bws <- as.data.frame(model, pars = c("theta1", "w"))
    coef <- mean(df_bws$theta1)
    coef_lb <- quantile(df_bws$theta1, 0.025)
    coef_ub <- quantile(df_bws$theta1, 0.975)
    
    # extract estimated weights
    weights <- colMeans(df_bws)[-1]
    weights_lb <- apply(df_bws, 2, quantile, probs = 0.025)[-1]
    weights_ub <- apply(df_bws, 2, quantile, probs = 0.975)[-1]
    p <- length(weights)
    
    # save results
    tmp_coef <- data.frame(method = method,
                           estimate = coef,
                           lb = coef_lb,
                           ub = coef_ub)
    
    tmp_weights <- data.frame(method = method,
                              mix_name = paste0("X",1:p),
                              estimate = weights,
                              lb = weights_lb,
                              ub = weights_ub)
    
    return(list(tmp_coef, tmp_weights))
  }
}
