### estimating beta ###
beta_first_mean <- apply(beta_first, 2, mean)
beta_final_mean <- apply(beta_final, 2, mean)

# BIAS
bias_first <- beta_first_mean - beta_actual
bias_final <- beta_final_mean - beta_actual

# SD
SD_first <- apply(beta_first, 2, sd)
SD_final <- apply(beta_final, 2, sd)

# SE
if (B != 0) {
  SE_first_mean <- apply(SE_first, 2, mean)
  SE_final_mean <- apply(SE_final, 2, mean)
}

# MSE
MSE_first <- bias_first^2 + SD_first^2 
MSE_final <- bias_final^2 + SD_final^2

# CP
if (B != 0) {
  # confidence interval of beta_first
  lower_first <- beta_first - 1.96 * SE_first # alpha = 0.05, C_alpha = 1.96
  upper_first <- beta_first + 1.96 * SE_first 
  confi_interval_first <- data.frame(matrix(nrow = N, ncol = 2))
  colnames(confi_interval_first) <- c("l", "u")
  if (is_ident) {
    for (j in 1 : (d+1)) {
      confi_interval_first[, (((j-1)*2+1) : (j*2))] <- data.frame(l = lower_first[ ,j], u = upper_first[ ,j])
    }
  } else {
    for (j in 1 : (K*d+1)) {
      confi_interval_first[, (((j-1)*2+1) : (j*2))] <- data.frame(l = lower_first[ ,j], u = upper_first[ ,j])
    }
  }
  # confidence interval of beta_final
  lower_final <- beta_final - 1.96 * SE_final # alpha = 0.05, C_alpha = 1.96
  upper_final <- beta_final + 1.96 * SE_final
  confi_interval_final <- data.frame(matrix(nrow = N, ncol = 2))
  colnames(confi_interval_final) <- c("l", "u")
  if (is_ident) {
    for (j in 1 : (d+1)) {
      confi_interval_final[, (((j-1)*2+1) : (j*2))] <- data.frame(l = lower_final[ ,j], u = upper_final[ ,j])
    }
  } else {
    for (j in 1 : (K*d+1)) {
      confi_interval_final[, (((j-1)*2+1) : (j*2))] <- data.frame(l = lower_final[ ,j], u = upper_final[ ,j])
    }
  }
}
if (B != 0) {
  # coverage probability (first)
  CP_first <- c()
  if (is_ident) {
    for (j in 1 : (d+1)) {
      CP_first_j <- sum(beta_actual[j] >= confi_interval_first[, (((j-1)*2+1) : (j*2))][1] &
                          beta_actual[j] <= confi_interval_first[, (((j-1)*2+1) : (j*2))][2]) / N
      CP_first <- c(CP_first, CP_first_j)
    }
  } else {
    for (j in 1 : (K*d+1)) {
      CP_first_j <- sum(beta_actual[j] >= confi_interval_first[, (((j-1)*2+1) : (j*2))][1] &
                          beta_actual[j] <= confi_interval_first[, (((j-1)*2+1) : (j*2))][2]) / N
      CP_first <- c(CP_first, CP_first_j)
    } 
  }
  # coverage probability (final)
  CP_final <- c()
  if (is_ident) {
    for (j in 1 : (d+1)) {
      CP_final_j <- sum(beta_actual[j] >= confi_interval_final[, (((j-1)*2+1) : (j*2))][1] &
                          beta_actual[j] <= confi_interval_final[, (((j-1)*2+1) : (j*2))][2]) / N
      CP_final <- c(CP_final, CP_final_j)
    }
  } else {
    for (j in 1 : (K*d+1)) {
      CP_final_j <- sum(beta_actual[j] >= confi_interval_final[, (((j-1)*2+1) : (j*2))][1] &
                          beta_actual[j] <= confi_interval_final[, (((j-1)*2+1) : (j*2))][2]) / N
      CP_final <- c(CP_final, CP_final_j)
    }
  }
}

# identical
maft_part <- data.frame(CR = censor_rate,
                        Method = c("FA-MAFT", "FA-MAFT", "MAFT", "MAFT"),
                        beta = c("beta1", "beta2", "beta1", "beta2"),
                        BIAS = c(bias_final[-1], bias_first[-1]), 
                        SD = c(SD_final[-1], SD_first[-1]),
                        SE = c(SE_final_mean[-1], SE_first_mean[-1]),
                        MSE = c(MSE_final[-1], MSE_first[-1]),
                        CP = c(CP_final[-1], CP_first[-1])
)
write.csv(maft_part, file = "maft_part.csv", row.names = TRUE)

### specific
## boxplot
library(ggplot2)
library(ggpubr)

# bias_final_scen10 <- bias_final[-1]
# bias_first_scen10 <- bias_first[-1]
# SD_final_scen10 <- SD_final[-1]
# SD_first_scen10 <- SD_first[-1]
# SE_final_mean_scen10 <- SE_final_mean[-1]
# SE_first_mean_scen10 <- SE_first_mean[-1]
# MSE_final_scen10 <- MSE_final[-1]
# MSE_first_scen10 <- MSE_first[-1]
# CP_final_scen10 <- CP_final[-1]
# CP_first_scen10 <- CP_first[-1]
# DSESD_final_scen10 <- (SE_final_mean-SD_final)[-1]
# DSESD_first_scen10 <- (SE_first_mean-SD_first)[-1]
# 
# 
# bias_final_scen20 <- bias_final[-1]
# bias_first_scen20 <- bias_first[-1]
# SD_final_scen20 <- SD_final[-1]
# SD_first_scen20 <- SD_first[-1]
# SE_final_mean_scen20 <- SE_final_mean[-1]
# SE_first_mean_scen20 <- SE_first_mean[-1]
# MSE_final_scen20 <- MSE_final[-1]
# MSE_first_scen20 <- MSE_first[-1]
# CP_final_scen20 <- CP_final[-1]
# CP_first_scen20 <- CP_first[-1]
# DSESD_final_scen20 <- (SE_final_mean-SD_final)[-1]
# DSESD_first_scen20 <- (SE_first_mean-SD_first)[-1]
# 
# bias_final_scen30 <- bias_final[-1]
# bias_first_scen30 <- bias_first[-1]
# SD_final_scen30 <- SD_final[-1]
# SD_first_scen30 <- SD_first[-1]
# SE_final_mean_scen30 <- SE_final_mean[-1]
# SE_first_mean_scen30 <- SE_first_mean[-1]
# MSE_final_scen30 <- MSE_final[-1]
# MSE_first_scen30 <- MSE_first[-1]
# CP_final_scen30 <- CP_final[-1]
# CP_first_scen30 <- CP_first[-1]
# DSESD_final_scen30 <- (SE_final_mean-SD_final)[-1]
# DSESD_first_scen30 <- (SE_first_mean-SD_first)[-1]

data_bias <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  Bias = c(bias_final_scen10, bias_first_scen10, bias_final_scen20, bias_first_scen20, bias_final_scen30, bias_first_scen30)
)
data_SD <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  SD = c(SD_final_scen10, SD_first_scen10, SD_final_scen20, SD_first_scen20, SD_final_scen30, SD_first_scen30)
)
data_SE <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  SE = c(SE_final_mean_scen10, SE_first_mean_scen10, SE_final_mean_scen20, SE_first_mean_scen20, SE_final_mean_scen30, SE_first_mean_scen30)
)
data_MSE <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  MSE = c(MSE_final_scen10, MSE_first_scen10, MSE_final_scen20, MSE_first_scen20, MSE_final_scen30, MSE_first_scen30)
)
data_abs_DCP095 <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  abs_DCP095 = abs(c(CP_final_scen10, CP_first_scen10, CP_final_scen20, CP_first_scen20, CP_final_scen30, CP_first_scen30) - 0.95)
)
data_abs_DSESD <- data.frame(
  CR = factor(rep(c("10%", "10%", "20%", "20%", "30%", "30%"), each=160)),
  Method  = factor(rep(c("FA-MAFT", "MAFT", "FA-MAFT", "MAFT", "FA-MAFT", "MAFT"), each=160)),
  abs_DSESD = abs(c(DSESD_final_scen10, DSESD_first_scen10, DSESD_final_scen20, DSESD_first_scen20, DSESD_final_scen30, DSESD_first_scen30))
)

# »æÖÆÏäÏßÍ¼
p1 <- ggplot(data_bias, aes(x=CR, y=Bias, fill = Method)) +
  geom_boxplot() +
  ylab("BIAS") +
  # ylim(-0.1, 0.1) + 
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

p2 <- ggplot(data_SD, aes(x=CR, y=SD, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

p3 <- ggplot(data_SE, aes(x=CR, y=SE, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

p4 <- ggplot(data_MSE, aes(x=CR, y=MSE, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

p5 <- ggplot(data_abs_DCP095, aes(x=CR, y=abs_DCP095, fill = Method)) +
  geom_boxplot() +
  ylab("|CP - 95%|") + 
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

p6 <- ggplot(data_abs_DSESD, aes(x=CR, y=abs_DSESD, fill = Method)) +
  geom_boxplot() +
  # ylab("|SD - SE|") +
  ylab(expression(paste("|", "SD - SE", "|"))) +
  scale_fill_manual(values = c("FA-MAFT" = "green", "MAFT" = "red")) +
  theme(legend.title = element_blank(), legend.position="none") +
  # theme(legend.title = element_blank(), legend.position = c(0.1, 0.95)) +
  geom_hline(yintercept = 0, linetype = "dashed")

ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, common.legend = T, legend = "bottom") 
