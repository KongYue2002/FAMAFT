library(aftgee)
library(Matrix) 
library(MASS)

library(foreach)
library(doParallel)
library(doRNG)  # control random numbers

source("functions.R")

##################################################################
#main procedure
##################################################################

### setting
##################################################################

## generate dat1
is_ident = FALSE
censor_rate = 0.2
DGP = "maft" # "famaft" or "maft" 
u_distr = "normal" # "logistic" or "normal"
Sigma_str = "lowrank" # "lowrank" or "AR"
n = 40
K = 80 # pay attention!!!
d = 2 # dimension of t(x_ki)
r = 3 # number of factor
if (is_ident) {
  beta = c(1, 1)
} else {
  beta = rep(c(1, 1, -1, 1, 1, -1, -1, -1), times = (K/4))
}
beta_actual = c(0, beta)

N = 100 # number of replications
B = 50 # times of resampling procedure 
B_var = 0 # If the aftgee function is not needed to compute beta_var, then B_var can be fixed to 0.

r_max = 5 #  used in the eigenvalue ratio method
iter_max = 15 # The maximum iteration count for arg_refine.
tol = 0.04

# Number of cores for parallel computing.
# n_cores <- detectCores() - 1  # Leave one core for the system.
n_cores <- 50 
##################################################################
### Store
##################################################################

censoring_rate <- c()
censoring_rate_margin <- data.frame(matrix(ncol = K))
colnames(censoring_rate_margin) <- paste("cen_rate", 1:K, sep = "")

### Store the original data.
T_N <- data.frame(matrix(nrow = K*n, ncol = N)) # The N columns correspond to the data T for I = 1, 2, ..., 100, respectively.
colnames(T_N) <- paste("T_rep", 1:N, sep = "")
Y_N <- data.frame(matrix(nrow = K*n, ncol = N)) # The N columns correspond to the data Y for I = 1, 2, ..., 100, respectively.
colnames(Y_N) <- paste("Y_rep", 1:N, sep = "")
if (is_ident) {
  X_df_N <- data.frame(matrix(ncol = d+1)) # The first K*n rows correspond to the data X_mat for I = 1.
  colnames(X_df_N) <- paste("x", 0:d, sep = "")
} else {
  X_df_N <- data.frame(matrix(ncol = K*d+1))
  colnames(X_df_N) <- paste("x", 0:(K*d), sep = "") # The first K*n rows correspond to the data X_mat for I = 1.
}

### Save the predicted data.
T_hat_final_N <- data.frame(matrix(nrow = K*n, ncol = N)) # The N columns correspond to the data T_hat_final for I = 1, 2, ..., 100, respectively.
colnames(T_hat_final_N) <- paste("T_hat_final_rep", 1:N, sep = "")
T_hat_first_N <- data.frame(matrix(nrow = K*n, ncol = N)) # The N columns correspond to the data T_hat_first for I = 1, 2, ..., 100, respectively.
colnames(T_hat_first_N) <- paste("T_hat_first_rep", 1:N, sep = "")
T_hat_TF_N <- data.frame(matrix(nrow = K*n, ncol = N)) # The N columns correspond to the data T_hat_TF for I = 1, 2, ..., 100, respectively.
colnames(T_hat_TF_N) <- paste("T_hat_TF_rep", 1:N, sep = "")
YT <- 0 # dim = c(K*n, 3). The 3 columns correspond Y_hat, T_0, and Y_hat - T_0, respectively.

# Store the regression parameter estimates; "first" denotes the results from the first iteration, and "final" denotes those from the last iteration.
if (is_ident) {
  beta_first <- beta_final <- SE_first <- SE_final <- data.frame(matrix(ncol = d+1))
  colnames(beta_first) <- colnames(beta_final) <- colnames(SE_first) <- colnames(SE_final) <- paste("beta", 0:d, sep = "")
} else {
  beta_first <- beta_final <- SE_first <- SE_final <- data.frame(matrix(ncol = K*d+1))
  colnames(beta_first) <- colnames(beta_final) <- colnames(SE_first) <- colnames(SE_final) <- paste("beta", 0:(K*d), sep = "")
}
ER_N <- data.frame(matrix(ncol = r_max)) # The eigenvalue ratios computed when estimating the number of factors in the last iteration.
colnames(ER_N) <- paste("ER_", 1:r_max, sep = "")
r_hat_1_histN <- list() 
F_first <- F_final <- data.frame(matrix(ncol = r_max)) # The number of columns is r_max.
colnames(F_first) <- colnames(F_final) <- paste("fac", 1:r_max, sep = "")
Lamb_first <- Lamb_final <- data.frame(matrix(ncol = r_max)) # The number of columns is r_max.
colnames(Lamb_first) <- colnames(Lamb_final) <- paste("load_on_fac", 1:r_max, sep = "")
common_error_first <- common_error_final <- c()
F_accur_first <- F_accur_final <- c()
Lamb_accur_first <- Lamb_accur_final <- c()

# Store all bootstrap results.
if (is_ident) {
  all_resam_beta_first <- all_resam_beta_final <- array(NA, dim = c(B, d+1, N))
} else {
  all_resam_beta_first <- all_resam_beta_final <- array(NA, dim = c(B, K*d+1, N))
}

# Store the conv from N runs of arg_refine.
convergence <- c()
convergence_step <- c()

### Do not consider the MAFT model; retain only the factor model.
r_hat_tf_N <- c()
common_error_tf_N <- c()
F_accur_tf_N <- c()
Lamb_accur_tf_N <- c()

#Store the computation times of the three methods.
time.famaft.N <- c()
time.maft.N <- c()
time.tf.N <- c()

##################################################################

cl <- makeCluster(n_cores)
registerDoParallel(cl)

registerDoRNG(123)

# Create an output directory.
dir.create("iteration_logs", showWarnings = FALSE)

results <- foreach (I = 1:N, .packages = c("Matrix", "MASS", "aftgee")) %dopar% {
  # Create a separate file for the current I.
  log_file <- sprintf("iteration_logs/I_%d.txt", I)
  # Optional: write the header.
  cat("I,J\n", file = log_file)
  
  dat1 <- datgen(n = n, K = K, r = r, beta = beta, is_ident = is_ident, censor_rate = censor_rate, DGP = DGP, u_distr = u_distr, Sigma_str = Sigma_str)
  while (exist_0col(dat1$data[c(-1, -2, -3)])) {
    dat1 <- datgen(n = n, K = K, r = r, beta = beta, is_ident = is_ident, censor_rate = censor_rate, DGP = DGP, u_distr = u_distr, Sigma_str = Sigma_str)
  }
  
  ### Save the original data.
  result_item <- list()
  result_item$censoring_rate <- dat1$censoring_rate
  result_item$censoring_rate_margin <- dat1$censoring_rate_margin
  result_item$T_data <- dat1$original$T_0
  result_item$Y_data <- dat1$original$Y_0
  
  if (is_ident) {
    result_item$X_data <- as.matrix(cbind(1, dat1$data[4 : (4+(d-1))]))
  } else {
    result_item$X_data <- as.matrix(cbind(1, dat1$data[4 : (4+(K*d-1))]))
  }
  
  ################################################################
  ### Do not consider the MAFT model; retain only the factor model.
  start.time.tf <- Sys.time()
  # Estimation 
  Y_mat <- matrix(dat1$original$Y_0, nrow = n, byrow = TRUE)
  f_num_tf <- f_num(Y_mat, r_max = r_max)
  r_hat_tf <- f_num_tf$r_hat
  fl_est_tf <- fl_est(Y_mat, r_hat = r_hat_tf)
  F_hat_tf <- fl_est_tf$F_hat
  Lamb_hat_tf <- fl_est_tf$Lamb_hat
  # 
  time.tf <- Sys.time() - start.time.tf
  result_item$time.tf <- as.numeric(time.tf, units = "secs")
  if (DGP == "famaft") {
    cfl_est_accur_tf <- cfl_est_accur(dat1$F_0, dat1$lambda, F_hat_tf, Lamb_hat_tf, r_hat = r_hat_tf)
  }
  # Store
  result_item$r_hat_tf <- r_hat_tf
  result_item$T_hat_TF <- c(t(as.matrix(F_hat_tf) %*% t(as.matrix(Lamb_hat_tf))))
  if (DGP == "famaft") {
    result_item$common_error_tf <- cfl_est_accur_tf$common_error
    result_item$F_accur_tf <- cfl_est_accur_tf$F_accur
    result_item$Lamb_accur_tf <- cfl_est_accur_tf$Lamb_accur
  }
  ################################################################
  start.time.famaft.1 <- Sys.time()
  refine_fit <- arg_refine(tol = tol, iter_max = iter_max, dat = dat1, r_max = r_max, is_ident = is_ident, B_var = B_var)
  time.famaft.1 <- Sys.time() - start.time.famaft.1
  time.famaft.1 <- as.numeric(time.famaft.1, units = "secs")
  
  result_item$YT <- refine_fit$YT
  result_item$time.maft <- refine_fit$time.maft
  result_item$time.famaft <- time.famaft.1 - refine_fit$time.store
  
  result_item$beta_first <- head(refine_fit$beta_hist, 1)
  ### fac_result
  result_item$F_first <- head(refine_fit$F_hist, n)
  result_item$Lamb_first <- head(refine_fit$Lamb_hist, K)
  if (DGP == "famaft") {
    result_item$common_error_first <- head(refine_fit$accuracy$common_error_hist, 1)
    result_item$F_accur_first <- head(refine_fit$accuracy$F_accur_hist, 1)
    result_item$Lamb_accur_first <- head(refine_fit$accuracy$Lamb_accur_hist, 1)
  }
  
  result_item$beta_final <- tail(refine_fit$beta_hist, 1) # tail(refine_fit$beta_hist, 1) Extract the last row of the data.frame.
  ### fac_result
  result_item$F_final <- tail(refine_fit$F_hist, n)
  result_item$Lamb_final <- tail(refine_fit$Lamb_hist, K)
  if (DGP == "famaft") {
    result_item$common_error_final <- tail(refine_fit$accuracy$common_error_hist, 1)
    result_item$F_accur_final <- tail(refine_fit$accuracy$F_accur_hist, 1) # tail(refine_fit$accuracy$F_accur_hist, 1) Extract the last element of the vector.
    result_item$Lamb_accur_final <- tail(refine_fit$accuracy$Lamb_accur_hist, 1)
  }
  
  result_item$r_hat_1_hist <- refine_fit$r_hat_1_hist
  result_item$ER <- refine_fit$ER
  
  ### Save the predicted data.
  r_hat_1_final <- refine_fit$r_hat_1_hist[length(refine_fit$r_hat_1_hist)]
  result_item$T_hat_final <- as.matrix(result_item$X_data) %*% as.numeric(result_item$beta_final) + c(t(as.matrix(result_item$F_final[, 1:r_hat_1_final]) %*% t(as.matrix(result_item$Lamb_final[, 1:r_hat_1_final]))))
  result_item$T_hat_first <- as.matrix(result_item$X_data) %*% as.numeric(result_item$beta_first)
  
  result_item$convergence <- refine_fit$conv
  result_item$convergence_step <- refine_fit$conv.step
  
  if (B != 0) {
    ### Compute SE_I via bootstrap.
    # Store the B beta estimates obtained from the B bootstrap resamples.
    if (is_ident) {
      resam_beta_first <- resam_beta_final <- data.frame(matrix(ncol = d+1))
      colnames(resam_beta_first) <- colnames(resam_beta_final) <- paste("beta", 0:d, sep = "")
    } else {
      resam_beta_first <- resam_beta_final <- data.frame(matrix(ncol = K*d+1))
      colnames(resam_beta_first) <- colnames(resam_beta_final) <- paste("beta", 0:(K*d), sep = "")
    }
    
    
    data <- dat1$data
    original <- dat1$original
    
    resam_beta <- for (J in 1:B) {
      # Record I and J
      cat(sprintf("%d,%d\n", I, J), file = log_file, append = TRUE)
      
      
      # Generate resampled sample
      re_index <- sample(1:n, replace = TRUE)
      while (length(unique(re_index)) <= 2) {
        re_index <- sample(1:n, replace = TRUE)
      }
      for (i in 1:n) {
        dat1$data[((i-1)*K+1) : (i*K), ] <- data[data$id == re_index[i], ]
        dat1$original[((i-1)*K+1) : (i*K), ] <- original[original$id == re_index[i], ]
      }
      while (exist_0col(dat1$data[c(-1, -2, -3)])) {
        re_index <- sample(1:n, replace = TRUE)
        while (length(unique(re_index)) <= 2) {
          re_index <- sample(1:n, replace = TRUE)
        }
        for (i in 1:n) {
          dat1$data[((i-1)*K+1) : (i*K),] <- data[data$id == re_index[i],]
          dat1$original[((i-1)*K+1) : (i*K), ] <- original[original$id == re_index[i], ]
        }
      }
      dat1$data$id <- rep(1:n, each = K)
      
      resample_fit <- arg_refine(tol = tol, iter_max = iter_max, dat = dat1, r_max = r_max, is_ident = is_ident, B_var = 0)
      
      resam_beta_first[J, ] <- as.numeric(head(resample_fit$beta_hist, 1))
      resam_beta_final[J, ] <- as.numeric(tail(resample_fit$beta_hist, 1))
    }
    result_item$resam_beta_first <- resam_beta_first
    result_item$resam_beta_final <- resam_beta_final
    
    result_item$SE_first <- apply(resam_beta_first, 2, sd) # 2 indicates columns
    result_item$SE_final <- apply(resam_beta_final, 2, sd) # 2 indicates columns
  }
  result_item 
}

# ============================================================================
# Combine the results.
# ============================================================================
for (I in 1:N) {
  censoring_rate[I] <- results[[I]]$censoring_rate
  censoring_rate_margin[I, ] <- results[[I]]$censoring_rate_margin
  T_N[, I] <- results[[I]]$T_data
  Y_N[, I] <- results[[I]]$Y_data
  
  if (is_ident) {
    X_df_N[(((I-1)*K*n+1) : (I*K*n)), ] <- results[[I]]$X_data
  } else {
    X_df_N[(((I-1)*K*n+1) : (I*K*n)), ] <- results[[I]]$X_data
  }
  
  time.tf.N[I] <- results[[I]]$time.tf
  T_hat_TF_N[, I] <- results[[I]]$T_hat_TF
  r_hat_tf_N[I] <- results[[I]]$r_hat_tf
  if (DGP == "famaft") {
    common_error_tf_N[I] <- results[[I]]$common_error_tf
    F_accur_tf_N[I] <- results[[I]]$F_accur_tf
    Lamb_accur_tf_N[I] <- results[[I]]$Lamb_accur_tf
  }
  
  YT <- YT + results[[I]]$YT
  
  time.maft.N[I] <- results[[I]]$time.maft
  time.famaft.N[I] <- results[[I]]$time.famaft
  
  beta_first[I, ] <- results[[I]]$beta_first
  beta_final[I, ] <- results[[I]]$beta_final
  
  F_first[(((I-1)*n+1) : (I*n)), ] <- results[[I]]$F_first
  F_final[(((I-1)*n+1) : (I*n)), ] <- results[[I]]$F_final
  Lamb_first[(((I-1)*K+1) : (I*K)), ] <- results[[I]]$Lamb_first
  Lamb_final[(((I-1)*K+1) : (I*K)), ] <- results[[I]]$Lamb_final
  
  if (DGP == "famaft") {
    common_error_first[I] <- results[[I]]$common_error_first
    common_error_final[I] <- results[[I]]$common_error_final
    F_accur_first[I] <- results[[I]]$F_accur_first
    F_accur_final[I] <- results[[I]]$F_accur_final
    Lamb_accur_first[I] <- results[[I]]$Lamb_accur_first
    Lamb_accur_final[I] <- results[[I]]$Lamb_accur_final
  }
  
  r_hat_1_histN[[I]] <- results[[I]]$r_hat_1_hist
  ER_N[I, ] <- results[[I]]$ER
  
  T_hat_final_N[, I] <- results[[I]]$T_hat_final
  T_hat_first_N[, I] <- results[[I]]$T_hat_first
  
  convergence[I] <- results[[I]]$convergence
  convergence_step[I] <- results[[I]]$convergence_step
  
  if (B != 0) {
    all_resam_beta_first[,, I] <- as.matrix(results[[I]]$resam_beta_first)
    all_resam_beta_final[,, I] <- as.matrix(results[[I]]$resam_beta_final)
    
    SE_first[I, ] <- results[[I]]$SE_first
    SE_final[I, ] <- results[[I]]$SE_final
  }
}
