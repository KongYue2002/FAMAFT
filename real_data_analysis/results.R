# p_value---------------------------
p_first <- 1 - apply(abs(beta_first/SE_first), 2, pnorm) + apply(-abs(beta_first/SE_first), 2, pnorm)
p_final <- 1 - apply(abs(beta_final/SE_final), 2, pnorm) + apply(-abs(beta_final/SE_final), 2, pnorm)

### beta_hat, p_value, and SE
beta <- data.frame(beta_first = t(beta_first), beta_final = t(beta_final))
p_SE <- data.frame(p_first, p_final, SE_first = t(SE_first), SE_final = t(SE_final))

### the estimated factors and loadings
F_hat_final <- F_final[,1:r_hat_1_final]
Lambda_hat_final <- Lamb_final[1:r_hat_1_final]

### the predicted failure time
T_hat_famaft <- T_hat_final_N[,1]
T_hat_maft <- T_hat_first_N[,1]
T_hat_tf <- T_hat_TF_N[,1]

Y_0 <- real_dat$original$Y_0
delta <- real_dat$data$status

### c_indices
library(survival)
fit_famaft <- concordance(Surv(Y_0, delta) ~ T_hat_famaft, reverse = FALSE)
fit_maft <- concordance(Surv(Y_0, delta) ~ T_hat_maft, reverse = FALSE)
fit_tf <- concordance(Surv(Y_0, delta) ~ T_hat_tf, reverse = FALSE)


c_indices <- data.frame(
  FA_MAFT = fit_famaft$concordance, 
  MAFT = fit_maft$concordance,
  TF = fit_tf$concordance
)

write.csv(c_indices, file = "c_indices.csv", row.names = TRUE)

### computation time ###
time_famaft <- mean(time.famaft.N)
time_maft <- mean(time.maft.N)
time_tf <- mean(time.tf.N)
time <- data.frame(FAMAFT = time_famaft, MAFT = time_maft, TF = time_tf) 

write.csv(time, file = "time.csv", row.names = TRUE)


