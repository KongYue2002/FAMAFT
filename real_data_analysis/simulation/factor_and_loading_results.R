### estimating common components, F and Lambda ###
if (DGP == "famaft") {
  AVE_CC <- c(mean(common_error_final), mean(common_error_tf_N))
  AVE_FS <- c(mean(F_accur_final), mean(F_accur_tf_N))
  AVE_FL <- c(mean(Lamb_accur_final), mean(Lamb_accur_tf_N))
  fac_part <- data.frame(AVE_CC = AVE_CC, AVE_FS = AVE_FS, AVE_FL = AVE_FL, row.names = c("FA-MAFT", "TF"))
}

write.csv(fac_part, file = "fac_part.csv", row.names = TRUE)

### selectting factors number r ###
r_hat_1_final_N <- c()
for (I in 1:N) {
  r_hat_1_final_N[I] <- r_hat_1_histN[[I]][length(r_hat_1_histN[[I]])]
}
propor <- c()
for (s in 1:r_max) {
  propor[s] <- sum(r_hat_1_final_N == s)/N
}
ER <- apply(ER_N, 2, mean)
fac_num <- data.frame(ER = ER, propor = propor)
rownames(fac_num) <- paste("r = ", 1:r_max, sep = "")

write.csv(t(fac_num), file = "fac_num.csv", row.names = TRUE)