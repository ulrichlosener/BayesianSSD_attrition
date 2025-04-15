#===============================================================================
# Function to calculate the effective sample size in a multilevel model with attrition
#===============================================================================

# Note that this function requires the packages "lme4" and "Matrix"

get_neff_attrition <- function(model, N, t.points, survival) {
  
  n <- length(t.points)
  
  # Extract random effects structure and residual variance
  reVar <- VarCorr(model)
  sigma2 <- sigma(model)^2
  
  # Extract subject-specific Z matrix
  Zt <- as.matrix(getME(model, "Zt")[1:2, 1:n])  # Transposed random-effects design matrix (Z')
  Z <- t(Zt)
  
  X_control <- matrix(data=c(rep(1, n), t.points, rep(0, n), rep(0, n)), nrow=n, ncol=4, dimnames=list(c(), c("intercept", "t", "treat", "t:treat")))
  X_treat <- matrix(data=c(rep(1, n), t.points, rep(1, n), t.points), nrow=n, ncol=4, dimnames=list(c(), c("intercept", "t", "treat", "t:treat")))
  
  # Construct D for one subject
  D <- as.matrix(bdiag(lapply(reVar, function(x) as(bdiag(x), "sparseMatrix"))))
  

  
  sum_mat <- sum_mat_indep <- matrix(rep(0, 16), nrow = 4, ncol = 4)
  
  # Backwards loop: from n (full data) to 1 (only first time point remains)
  # k represents the last time point of measurement
  for (k in n:1) {
    if (k == n) {
      # Full dataset (no attrition)
      Z_subset <- Z
      X_control_subset <- X_control
      X_treat_subset <- X_treat
    } else {
      # Subset to keep only rows 1:k
      Z_subset <- Z[1:k, , drop = FALSE]
      X_control_subset <- X_control[1:k, , drop = FALSE]
      X_treat_subset <- X_treat[1:k, , drop = FALSE]
    }
    # Construct V for this specific attrition pattern
    V_subset <- Z_subset %*% D %*% t(Z_subset) + sigma2 * Diagonal(nrow(Z_subset))
    W_subset <- Diagonal(x = diag(V_subset))
    
    # Take inverses of V and W
    V_subset_inv <- solve(V_subset)
    W_subset_inv <- solve(W_subset)
    
    # Compute sum_mat and sum_mat_indep weighted by S(t)
    sum_mat_control <- t(X_control_subset) %*% V_subset_inv %*% X_control_subset
    sum_mat_treat <- t(X_treat_subset) %*% V_subset_inv %*% X_treat_subset
    sum_mat <- sum_mat + (survival[k] * as.matrix(sum_mat_control + sum_mat_treat))
    
    sum_mat_control_indep <- t(X_control_subset) %*% W_subset_inv %*% X_control_subset
    sum_mat_treat_indep <- t(X_treat_subset) %*% W_subset_inv %*% X_treat_subset
    sum_mat_indep <- sum_mat_indep + (survival[k] * as.matrix(sum_mat_control_indep + sum_mat_treat_indep))
    
  }
  
  # if sum_mat is near-singular use pseudoinverse
  if (rcond(sum_mat) < 1e-10) { 
    var_beta_hat <- MASS::ginv(sum_mat)
  } else {
    var_beta_hat <- solve(sum_mat)
  }
  
  if (rcond(sum_mat) < 1e-10) { 
    var_betahat_indep <- MASS::ginv(sum_mat_indep)
  } else {
    var_betahat_indep <- solve(sum_mat_indep)
  }
  
  # Compute weight w and effective sample size N_eff
  w <- var_betahat_indep / var_beta_hat
  N_eff <- w * N * n
  
  return(N_eff[4, 4])
}


