### BWS ###
library(bws)
library(rstan)

### suppose you have outcome y, exposures x, and covariates z
fit.bws <- bws(iter = 2000, y = y, X = x, Z = z, family = "gaussian", 
               show_messages = FALSE)
# display theta1 and w
print(fit.bws, pars = c("w", "theta1"))
# visualization
rstan::plot(fit.bws, pars = c("w", "theta1"))

# extract estimated mixture effect with CI
df_bws <- as.data.frame(fit.bws, pars = c("theta1", "w"))
coef <- mean(df_bws$theta1)
coef_lb <- quantile(df_bws$theta1, 0.025)
coef_ub <- quantile(df_bws$theta1, 0.975)

# extract estimated weights with CIs
weights <- colMeans(df_bws)[-1]
weights_lb <- apply(df_bws, 2, quantile, probs = 0.025)[-1]
weights_ub <- apply(df_bws, 2, quantile, probs = 0.975)[-1]
