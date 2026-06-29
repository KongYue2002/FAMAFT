##################################################################
#data generation
##################################################################

datgen <- function(n, K , r, beta, is_ident, censor_rate, DGP, u_distr, Sigma_str) {#n subjects, K events, r factors
  x1_data <- rbinom(K * n, size = 1, prob = 0.5) # the 1st covariate X_1ki
  x2_data <- rnorm(K * n, mean = 0, sd = 1) # X_2ki; dimension of covariates d = 2 
  X_data <- cbind(x1_data, x2_data)
  if (is_ident) {
    X <- cbind(x1_data, x2_data)
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
  
  if (DGP == "famaft") {
    ## factor model
    lamb <- matrix(runif(K * r, min = -1, max = 1), nrow = K, ncol = r)
    e <- c()
    F_0 <- matrix(NA, nrow = n, ncol = r)
    for (i in 1:n) {
      f_i <- rnorm(r, mean = 0, sd = 1)
      if (u_distr == "logistic") {
        u_i <- rlogis(K, location = 0, scale = 0.1)
      }
      if (u_distr == "normal") {
        u_i <- rnorm(K, mean = 0, sd = 0.1)
      }
      e[((i - 1) * K + 1) : (i * K)] <- lamb %*% f_i + u_i # e_i <- lamb * f_i + u_i # %*% Matrix inner product
      F_0[i,] <- f_i
    }
  }
  if (DGP == "maft") {
    if (Sigma_str == "lowrank") {
      A <- matrix(0, nrow = K, ncol = r)
      A[,1] <- 0.8
      A[,2] <- 0.6
      A[,3] <- 0.4
      Sigma <- A %*% t(A)
    }
    if (Sigma_str == "AR") {
      rho <- 0.85
      Sigma <- rho^abs(outer(1:K, 1:K, "-"))
    }
    U <- mvrnorm(n = n, mu = rep(0, K), Sigma = Sigma) # Each row of U is a K-dimensional normal random vector. Generate n multivariate normal random numbers.
    u <- c(t(U))
    e <- u
  }
  ## logT_i = X_i * beta + epsilon
  tt <- exp(X %*% beta + e) 
  
  if (is_ident) {
    ### apply censoring to log(tt)
    if (censor_rate == 0.1) {
      log_cen <- runif(K*n, 0, 8.8) # censoring_rate = 10%
    }
    if (censor_rate == 0.2) {
      log_cen <- runif(K*n, 0, 4.4) # censoring_rate = 20%
    }
    if (censor_rate == 0.3) {
      log_cen <- runif(K*n, 0, 2.8) # censoring_rate = 30%
    }
    delta <- log(tt) <= log_cen
    cen_rate <- length(delta[delta == FALSE]) / (K*n)
    del_mat <- matrix(delta, nrow = n, ncol = K, byrow = TRUE)
    cen_rate_margin <- 1 - apply(del_mat, 2, mean) # dim=c(K£¬1), 2 in apply() indicates columns
    C_0 <- log_cen
  } else {
    ### apply censoring to each of the K margins of log(tt)
    if (censor_rate == 0.1) {
      c <- rep(c(8.8, 3.65, 8.8, 3.65), times = (K/4)) # censoring_rate = 10%
    }
    if (censor_rate == 0.2) {
      c <- rep(c(4.4, 1.6, 4.4, 1.6), times = (K/4)) # censoring_rate = 20%
    }
    if (censor_rate == 0.3) {
      c <- rep(c(2.8, 0.547, 2.8, 0.547), times = (K/4)) # censoring_rate = 30%  
    }
    log_cen_mat <- matrix(nrow = n, ncol = K)
    T_mat <- matrix(log(tt), nrow = n, byrow = TRUE)
    del_mat <- matrix(nrow = n, ncol = K)
    for (k in 1:K) {
      # log_cen_mat[k] = (C_k1,..., C_kn), k = 1,...,K
      log_cen_mat[,k] <- runif(n, 0, c[k]) # dim = c(n, 1)
      del_mat[,k] <- T_mat[,k] <= log_cen_mat[,k]
    }
    cen_rate_margin <- 1 - apply(del_mat, 2, mean) # dim=c(K£¬1), 2 in apply() indicates columns
    cen_rate <- mean(cen_rate_margin)
    C_0 = c(t(log_cen_mat)) # dim = c(K*n, 1)  
  }
  
  
  da_li <- list(
    original = data.frame(T_0 = log(tt),
                          C_0 = C_0,
                          Y_0 = pmin(log(tt), C_0),
                          epsilon_0 = e,
                          id = rep(1:n, each = K)),
    censoring_rate_margin = cen_rate_margin,
    censoring_rate = cen_rate,
    covariates = X,
    data = data.frame(Time = pmin(tt, exp(C_0)),
                      status = 1 * (log(tt) <= C_0),
                      id = rep(1:n, each = K)))
  if (is_ident) {
    da_li$data[4 : (4+(d-1))] <- X
  } else {
    da_li$data[4 : (4+(K*d-1))] <- X # X
  }
  if (DGP == "famaft") {
    da_li$F_0 <- F_0
    da_li$lambda <- lamb
  }
  da_li
}

### define a function that checks whether df has any zero columns, and returns TRUE if so
exist_0col <- function (df) {
  is_0col <- c()
  for (j in 1:ncol(df)) {
    is_0col[j] <- all(df[j] == 0)
  }
  any(is_0col)
}

### estimate factors' number by eigenvalue ratio method
f_num <- function (data_mat, r_max) {
  cov <- cov(data_mat)
  eig <- eigen(cov)
  eig_va <- eig$values
  eig_vec <- eig$vectors
  
  eig_ratio <- c()
  for (s in 1:r_max) {
    eig_ratio[s] <- eig_va[s]/eig_va[s+1]
  }
  r_hat <- which.max(eig_ratio) # use the which.max function to obtain the position of the maximum value in a vector
  list(r_hat = r_hat, eig_ratio = eig_ratio)
}


fl_est <- function (data_mat, r_hat) {
  ### estimate factors and factor loadings by least squares (equivalent to PCA)
  col_mean <- apply(data_mat, MARGIN = 2, mean)
  data_cen <- t(t(data_mat) - col_mean) ## remove the mean of each variable of data_mat
  F_hat <- sqrt(n) * eigen(data_cen %*% t(data_cen))$vector[, 1:r_hat, drop = FALSE] # dim = c(n, r)
  Lamb_hat <- t(data_cen) %*% F_hat / n # dim = c(K, r)
  
  list(F_hat = F_hat, Lamb_hat = Lamb_hat)
}

### the estimation error for common components, loading matrices and factor score matrices
cfl_est_accur <-  function (F_0, lambda, F_hat, Lamb_hat, r_hat) {
  # Gram¨CSchmidt orthogonalization
  # QR decomposition
  qr_F_0 <- qr(F_0)
  qr_lambda <- qr(lambda)
  qr_F_hat <- qr(F_hat)
  qr_Lamb_hat <- qr(Lamb_hat)
  # extract the Q matrix from the QR decomposition
  F_0_orth <- qr.Q(qr_F_0)
  lambda_orth <- qr.Q(qr_lambda)
  F_hat_orth <- qr.Q(qr_F_hat)
  Lamb_hat_orth <- qr.Q(qr_Lamb_hat)
  
  common_error <- sum((Lamb_hat %*% t(F_hat) - lambda %*% t(F_0))^2) / sum((lambda %*% t(F_0))^2)
  F_accur <- sqrt(1 - sum(diag(F_hat_orth %*% t(F_hat_orth) %*%  F_0_orth %*% t( F_0_orth))) / r_hat)
  Lamb_accur <- sqrt(1 - sum(diag(Lamb_hat_orth %*% t(Lamb_hat_orth) %*%  lambda_orth %*% t( lambda_orth))) / r_hat)
  
  list(common_error = common_error, F_accur = F_accur, Lamb_accur = Lamb_accur)
}

# Compute the LS integral; f is the function F(x); 
# a and b are the endpoints of the integration interval; 
# n is the number of subintervals into which the integration interval is divided.
LS_integrate <- function (f, a, b, n = 50) { 
  integral_value <- 0
  dx <- (b - a) / n
  points <- seq(a, b, dx) # points[i] is x_{i-1}; seq(0, 1, 0.1) = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
  for (i in 1:n) {
    integral_value <-  integral_value + (points[i] + points[i+1])/2 * (f(points[i+1]) - f(points[i]))
  }
  integral_value
}

###############################################
### Compute the IPCW estimate of the observed failure time Y. 
# (The LS_integrate function must be defined before the IPCW function.)
###############################################
IPCW <- function (Y, X_mat, beta_hat, delta) {
  # Y: a Kn-dimensional vector; X_mat: a matrix with Kn rows and either Kd+1 or d+1 columns; beta_hat: a vector of dimension Kd+1 or d+1; delta_mat: an (n,K) matrix.
  
  e_hat <- Y - X_mat%*%beta_hat 
  e_hat_mat <- matrix(e_hat, nrow = n, byrow = TRUE)
  
  delta_mat <- matrix(delta, nrow = n, byrow = TRUE)
  
  ### Estimation of the distribution function of epsilon_ki, k = 1, ..., K
  ## Define a vector containing the strings of function definitions.
  Fun_hat_strvec <- c() # Define a vector to store strings.
  # Use a for loop to create a character vector.
  for (k in 1:K) {
    # Create a string that contains the value of the loop variable.
    Fun_hat_str <- paste0(
      "function (t) {
      su <- 1
      for (i in 1:n) {
        if (e_hat_mat[i,", k, "] < t) {
          su <- su * (1 - delta_mat[i,", k, "] / sum(e_hat_mat[,", k, "] >= e_hat_mat[i,", k, "]))
        }
      }
      1 - su
    }"
    )
    Fun_hat_strvec <- c(Fun_hat_strvec, Fun_hat_str) # Append the newly created string to the string vector.
  }
  ## Parse the string and define the function.
  Fun_hat <- lapply(Fun_hat_strvec, function(s) {
    eval(parse(text = s))
  })
  
  e_ipcw_mat <- matrix(data = NA, nrow = n, ncol = K)
  for (k in 1:K) {
    for (i in 1:n) {
      Ee <- LS_integrate(f = Fun_hat[[k]], a = e_hat_mat[i,k], b = max(e_hat_mat[,k]) + 0.1, n = 50) # Adding 0.1 prevents a and b from being equal.
      e_ipcw_mat[i,k] <- Ee/(1 - Fun_hat[[k]](e_hat_mat[i,k]))
    }
  }
  
  e_ipcw <- c(t(e_ipcw_mat))
  Y_hat <- delta * Y + (1 - delta) * (e_ipcw + X_mat %*% beta_hat)
  
  Y_hat
}




##################################################################
#iteration steps to update beta, f and lambda until convergence
##################################################################

arg_refine <- function (tol , iter_max, dat, r_max, is_ident, B_var) {
  M <- 1
  conv <- 0
  
  start.time.maft <- Sys.time()
  #######################################
  if (is_ident) {
    x_fm <- paste("V", 4:(4+(d-1)), collapse = " + ")
  } else {
    # Create a string containing all the covariates.
    x_fm <- paste("V", 4:(4+(K*d-1)), collapse = " + ")
  }
  # Create a formula string.
  xy_fm <- paste("Surv(Time, status) ~", x_fm)
  # Use gsub() to remove all spaces.
  xy_fm <- gsub("\\s", "", xy_fm)
  # Convert the formula string to a formula object.
  fm <- as.formula(xy_fm)
  
  fit <- aftgee(fm, data = dat$data, id = id, corstr = "ind", B = B_var)
  beta_hat <- fit$coef.res  # a vector of point estimates
  #######################################
  time.maft <- Sys.time() - start.time.maft
  
  delta <- dat$data$status
  Y <- dat$original$Y_0
  if (is_ident) {
    X_mat <- as.matrix(cbind(1, dat$data[4 : (4+(d-1))]))
  } else {
    X_mat <- as.matrix(cbind(1, dat$data[4 : (4+(K*d-1))]))
  }
  
  Y_hat <- IPCW(Y = Y, X_mat = X_mat, beta_hat = beta_hat, delta = delta)
  epsilon_hat <- Y_hat - X_mat %*% beta_hat
  epsilon_hat_mat <- matrix(epsilon_hat, nrow = n, byrow = TRUE)
  
  # Estimate the factor number.
  f_num_famaft <- f_num(epsilon_hat_mat, r_max = r_max)
  r_hat_1 <- f_num_famaft$r_hat
  eig_ratio_1 <- f_num_famaft$eig_ratio
  
  ### estimate factors and factor loadings
  epsilon_hat_fl_est <- fl_est(epsilon_hat_mat, r_hat = r_hat_1)
  F_hat_IPCW <- epsilon_hat_fl_est$F_hat # dim = c(n, r)
  Lamb_hat_IPCW <- epsilon_hat_fl_est$Lamb_hat # dim = c(K, r)
  
  start_time.store.1 <- Sys.time()
  ############################################################
  # Store the parameter estimates.
  if (is_ident) {
    beta_hist <- data.frame(matrix(ncol = d+1))
    colnames(beta_hist) <- paste("beta", 0:d, sep = "")
    beta_var <- data.frame(matrix(ncol = d+1))
    colnames(beta_var) <- paste("beta", 0:d, sep = "")
  } else {
    beta_hist <- data.frame(matrix(ncol = K*d+1))
    colnames(beta_hist) <- paste("beta", 0:(K*d), sep = "")
    beta_var <- data.frame(matrix(ncol = K*d+1))
    colnames(beta_var) <- paste("beta", 0:(K*d), sep = "") 
  }
  r_hat_1_hist <- c()
  F_hist <- data.frame(matrix(ncol = r_max)) # The number of columns is r_max.
  colnames(F_hist) <- paste("fac", 1:r_max, sep = "")
  Lamb_hist <- data.frame(matrix(ncol = r_max)) # The number of columns is r_max.
  colnames(Lamb_hist) <- paste("load_on_fac", 1:r_max, sep = "")
  
  beta_hist[M,] <- fit$coef.res
  if (B_var != 0) {
    beta_var[1,] <- diag(fit$var.res)
  }
  r_hat_1_hist[M] <- r_hat_1
  F_hist[(((M-1)*n+1) : (M*n)), 1:r_hat_1] <- F_hat_IPCW
  Lamb_hist[(((M-1)*K+1) : (M*K)), 1:r_hat_1] <- Lamb_hat_IPCW
  
  # Store the estimation accuracy of the factor model parameters.
  common_error_hist <- c() # normalized estimation error for common components in terms of matrix Frobenius norm
  F_accur_hist <- c()
  Lamb_accur_hist <- c()
  if (DGP == "famaft") {
    epsilon_hat_cfl_est_accur <- cfl_est_accur(dat$F_0, dat$lambda, F_hat_IPCW, Lamb_hat_IPCW, r_hat = r_hat_1)
    common_error_hist[M] <- epsilon_hat_cfl_est_accur$common_error
    F_accur_hist[M] <- epsilon_hat_cfl_est_accur$F_accur
    Lamb_accur_hist[M] <- epsilon_hat_cfl_est_accur$Lamb_accur
  }
  ############################################################
  time.store.1 <- Sys.time() - start_time.store.1
  
  Y_0 <- dat$original$Y_0 # Store the initial (censored log-transformed) times, i.e., Y_i, i = 1, 2, ..., n (not involved in the loop).
  
  while (M <= iter_max) {
    for (i in 1:n) {
      dat$data$Time[((i - 1)*K + 1) : (i*K)] <-  exp(Y_0[((i - 1)*K + 1) : (i*K)] - Lamb_hat_IPCW %*% F_hat_IPCW[i, ])
      # log(Time0[((i - 1)*K + 1) : (i*K)]) is Y_i
      # Lamb_hat_IPCW is Lambda^hat_(M)
      # F_hat_IPCW is F^hat_(M)
    }
    
    fit <- aftgee(fm, data = dat$data, id = id, corstr = "ind", binit =as.numeric(beta_hist[M,]), B = 0)
    beta_hat <- fit$coef.res  # a vector of point estimates
    
    Y_hat <- IPCW(Y = Y, X_mat = X_mat, beta_hat = beta_hat, delta = delta)
    epsilon_hat <- Y_hat - X_mat %*% beta_hat
    epsilon_hat_mat <- matrix(epsilon_hat, nrow = n, byrow = TRUE)
    
    # Estimate the factor number.
    f_num_famaft <- f_num(epsilon_hat_mat, r_max = r_max)
    r_hat_1 <- f_num_famaft$r_hat
    eig_ratio_1 <- f_num_famaft$eig_ratio
    
    ### estimate factors and factor loadings
    epsilon_hat_fl_est <- fl_est(epsilon_hat_mat, r_hat = r_hat_1)
    F_hat_IPCW <- epsilon_hat_fl_est$F_hat # dim = c(n, r)
    Lamb_hat_IPCW <- epsilon_hat_fl_est$Lamb_hat # dim = c(K, r)
    
    start_time.store.2 <- Sys.time()
    ############################################################
    # Store the parameter estimates.
    beta_hist[M+1,] <- fit$coef.res # a vector of point estimates
    r_hat_1_hist[M+1] <- r_hat_1
    F_hist[((((M+1)-1)*n+1) : ((M+1)*n)), 1:r_hat_1] <- F_hat_IPCW
    Lamb_hist[((((M+1)-1)*K+1) : ((M+1)*K)), 1:r_hat_1] <- Lamb_hat_IPCW
    
    # Store the estimation accuracy of the factor model parameters.
    if (DGP == "famaft") {
      epsilon_hat_cfl_est_accur <- cfl_est_accur(dat$F_0, dat$lambda, F_hat_IPCW, Lamb_hat_IPCW, r_hat = r_hat_1)
      common_error_hist[M+1] <- epsilon_hat_cfl_est_accur$common_error
      F_accur_hist[M+1] <- epsilon_hat_cfl_est_accur$F_accur
      Lamb_accur_hist[M+1] <- epsilon_hat_cfl_est_accur$Lamb_accur
    }
    ############################################################
    time.store.2 <- Sys.time() - start_time.store.2
    
    # Convergence criterion: the distance between two consecutive beta_hat estimates is less than tol.
    # norm <- sqrt(as.matrix((beta_hist[M,] - beta_hist[M+1,])) %*% t(as.matrix((beta_hist[M,] - beta_hist[M+1,]))))
    norm <- max(abs(beta_hist[M,] - beta_hist[M+1,]))
    if (norm < tol) {
      conv <- 1
      if (B_var != 0) {
        fit <- aftgee(fm, data = dat$data, id = id, corstr = "ind", binit =as.numeric(beta_hist[M,]), B = B_var)
        beta_var[2,] <- diag(fit$var.res)
      }
      break
    }
    M <- M + 1
  }
  if (conv == 1) {
    conv.step <- M + 1
  } else {
    conv.step <- M
  }
  refine_result <- list(beta_hist = beta_hist,
                        beta_var = beta_var,
                        F_hist = F_hist, 
                        Lamb_hist = Lamb_hist,
                        r_hat_1_hist = r_hat_1_hist,
                        ER = eig_ratio_1, # The eigenvalue ratios computed when estimating the number of factors in the last iteration.
                        accuracy = data.frame(common_error_hist = common_error_hist,
                                              F_accur_hist = F_accur_hist,
                                              Lamb_accur_hist = Lamb_accur_hist),
                        YT = cbind(Y_hat, dat$original$T_0, Y_hat - dat$original$T_0),
                        time.maft = as.numeric(time.maft, units = "secs"),
                        time.store = as.numeric(time.store.1 + time.store.2, units = "secs"),
                        conv = conv, # Indicates whether convergence is achieved; conv = 1 means convergence.
                        conv.step = conv.step ) # Denotes the iteration count.
  refine_result
}


