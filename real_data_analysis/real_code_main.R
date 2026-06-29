library(aftgee)
library(Matrix)

library(foreach)
library(doParallel)
library(doRNG)  # control random numbers

source("functions.R")


#########################################
#real data################################################
#########################################


###setting###############################
is_ident = FALSE

DGP = "maft" # "famaft" or "maft" 

n <- 30
K <- 76
d <- 4

N = 1
B = 0 # times of resampling procedure 
B_var = 0

r_max = 5 # used in the eigenvalue ratio method
iter_max = 15 # The maximum iteration count for arg_refine.
tol = 0.04

# Number of cores for parallel computing.
# n_cores <- detectCores() - 1  # Leave one core for the system.


exp_Y0_init <- read.csv("data\\exp_Y0.csv")$exp_Y0
delta_init <- read.csv("data\\delta.csv")$delta
X_dat <- read.csv("data\\X.csv",header = TRUE)

exp_Y0_mat <- matrix(exp_Y0_init, nrow = n, byrow = FALSE)
delta_mat <- matrix(delta_init, nrow = n, byrow = FALSE)
X_data <- matrix(nrow = n*K, ncol = d)
X_data[,1] <- X_dat[,1] # No perturbation needs to be added to x1.
X_data[,2:d] <- as.matrix(X_dat[1:n,2:d])[rep(1:n, each = K), ]
X_data[,2:d] <- X_data[,2:d] + matrix(rnorm(K*n*(d-1), mean = 0, sd = 0.01), nrow = K*n, ncol = d-1)
if (is_ident) {
  X <- X_data
} else {
  X <- matrix(nrow = K*n, ncol = K*d)
  for (i in 1:n) {
    X_data_i <- X_data[((i - 1) * K + 1) : (i * K),]
    X_i <- matrix()
    for ( k in 1:K) {
      x_ki <- matrix(X_data_i[k,], nrow = 1, byrow = TRUE)  
      X_i <- as.matrix(bdiag(X_i, x_ki))
    }
    X[((i - 1) * K + 1) : (i * K), ] <- as.matrix(X_i[-1, -1]) # pay attention!!! K can't be 1
  }
}



exp_Y0 <- c(t(exp_Y0_mat))
delta <- c(t(delta_mat))            

cen_rate_margin <- 1 - apply(delta_mat, 2, mean)
cen_rate = 1 - mean(delta)

realdata <- data.frame(Time = exp_Y0, 
                       status = delta, 
                       id = rep(1:n, each = K))
if (is_ident) {
  realdata[4 : (4+(d-1))] <- X
} else {
  realdata[4 : (4+(K*d-1))] <- X 
}

real_dat <- list(F_0 = NA,
                 lambda = NA,
                 original = data.frame(T_0 = NA,
                                       C_0 = NA,
                                       Y_0 = log(exp_Y0),
                                       epsilon_0 = NA,
                                       id = rep(1:n, each = K)),
                 censoring_rate_margin = cen_rate_margin,
                 censoring_rate = cen_rate,
                 covariates = X,
                 data = realdata
)

#########################################
#real main################################################
#########################################

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
YT <- 0 # dim = c(K*n, 3) The 3 columns correspond Y_hat, T_0, and Y_hat - T_0, respectively.

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


for (I in 1:N) {
  print(I)
  set.seed(I * 123)
  
  censoring_rate <- c(censoring_rate, real_dat$censoring_rate)
  censoring_rate_margin[I,] <- real_dat$censoring_rate_margin
  
  ### Save the original data.
  Y_N[,I] <- real_dat$original$Y_0
  if (is_ident) {
    X_df_N[(((I-1)*K*n+1) : (I*K*n)), ] <- as.matrix(cbind(1, real_dat$data[4 : (4+(d-1))]))
  } else {
    X_df_N[(((I-1)*K*n+1) : (I*K*n)), ] <- as.matrix(cbind(1, real_dat$data[4 : (4+(K*d-1))]))
  }
  
  ################################################################
  ### Do not consider the MAFT model; retain only the factor model.
  start.time.tf <- Sys.time()
  # Estimation 
  Y_mat <- matrix(real_dat$original$Y_0, nrow = n, byrow = TRUE)
  f_num_tf <- f_num(Y_mat, r_max = r_max)
  r_hat_tf <- f_num_tf$r_hat
  fl_est_tf <- fl_est(Y_mat, r_hat = r_hat_tf)
  F_hat_tf <- fl_est_tf$F_hat
  Lamb_hat_tf <- fl_est_tf$Lamb_hat
  # 
  time.tf <- Sys.time() - start.time.tf
  time.tf.N[I] <- as.numeric(time.tf, units = "secs")
  # Store
  r_hat_tf_N <- c(r_hat_tf_N, r_hat_tf)
  T_hat_TF_N[,I] <- c(t(as.matrix(F_hat_tf) %*% t(as.matrix(Lamb_hat_tf))))
  ################################################################
  start.time.famaft.1 <- Sys.time()
  refine_fit <- arg_refine(tol = tol, iter_max = iter_max, dat = real_dat, r_max = r_max, is_ident = is_ident, B_var = B_var)
  time.famaft.1 <- Sys.time() - start.time.famaft.1
  time.famaft.1 <- as.numeric(time.famaft.1, units = "secs")
  
  YT <- YT + refine_fit$YT
  time.maft.N[I] <- refine_fit$time.maft
  time.famaft.N[I] <- time.famaft.1 - refine_fit$time.store
  
  beta_first[I, ] <- head(refine_fit$beta_hist, 1)
  ### fac_result
  F_first[(((I-1)*n+1) : (I*n)), ] <- head(refine_fit$F_hist, n)
  Lamb_first[(((I-1)*K+1) : (I*K)), ] <- head(refine_fit$Lamb_hist, K)
  
  beta_final[I, ] <- tail(refine_fit$beta_hist, 1) # tail(refine_fit$beta_hist, 1) Extract the last row of the data.frame.
  ### fac_result
  F_final[(((I-1)*n+1) : (I*n)), ] <- tail(refine_fit$F_hist, n)
  Lamb_final[(((I-1)*K+1) : (I*K)), ] <- tail(refine_fit$Lamb_hist, K)
  
  r_hat_1_histN <- c(r_hat_1_histN, list(refine_fit$r_hat_1_hist))
  ER_N[I,] <- refine_fit$ER
  
  ### Save the predicted data.
  r_hat_1_final <- refine_fit$r_hat_1_hist[length(refine_fit$r_hat_1_hist)]
  T_hat_final_N[,I] <- as.matrix(X_df_N[(((I-1)*K*n+1) : (I*K*n)), ]) %*% as.numeric(beta_final[I, ]) + c(t(as.matrix(F_final[(((I-1)*n+1) : (I*n)), 1:r_hat_1_final]) %*% t(as.matrix(Lamb_final[(((I-1)*K+1) : (I*K)), 1:r_hat_1_final]))))
  T_hat_first_N[,I] <- as.matrix(X_df_N[(((I-1)*K*n+1) : (I*K*n)), ]) %*% as.numeric(beta_first[I, ])
  
  convergence <- c(convergence, refine_fit$conv)
  convergence_step <- c(convergence_step, refine_fit$conv.step)
  
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
    
    
    data <- real_dat$data
    original <- real_dat$original
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    registerDoRNG(I * 456)
    
    resam_beta <- foreach(J=1:B, .combine='rbind', .packages=c("aftgee")) %dopar% {
      # Generate resampled sample
      re_index <- sample(1:n, replace = TRUE)
      while (length(unique(re_index)) <= 2) {
        re_index <- sample(1:n, replace = TRUE)
      }
      for (i in 1:n) {
        real_dat$data[((i-1)*K+1) : (i*K), ] <- data[data$id == re_index[i], ]
        real_dat$original[((i-1)*K+1) : (i*K), ] <- original[original$id == re_index[i], ]
      }
      while (exist_0col(real_dat$data[c(-1, -2, -3)])) {
        re_index <- sample(1:n, replace = TRUE)
        while (length(unique(re_index)) <= 2) {
          re_index <- sample(1:n, replace = TRUE)
        }
        for (i in 1:n) {
          real_dat$data[((i-1)*K+1) : (i*K),] <- data[data$id == re_index[i],]
          real_dat$original[((i-1)*K+1) : (i*K), ] <- original[original$id == re_index[i], ]
        }
      }
      real_dat$data$id <- rep(1:n, each = K)
      
      resample_fit <- arg_refine(tol = tol, iter_max = iter_max, dat = real_dat, r_max = r_max, is_ident = is_ident, B_var = 0)
      
      return(as.numeric(c(head(resample_fit$beta_hist, 1), tail(resample_fit$beta_hist, 1))))
    }
    stopCluster(cl)
    
    if (is_ident) {
      resam_beta_first <- resam_beta[,1:(d+1)]
      resam_beta_final <- resam_beta[,(d+2):ncol(resam_beta)]
    } else {
      resam_beta_first <- resam_beta[,1:(K*d+1)]
      resam_beta_final <- resam_beta[,(K*d+2):ncol(resam_beta)]
    }
    SE_first[I, ] <- apply(resam_beta_first, 2, sd) # 2 indicates columns
    SE_final[I, ] <- apply(resam_beta_final, 2, sd) # 2 indicates columns
  }
}