#===============================================================================
# Function to calculate the effective sample size in a multilevel model with attrition
#===============================================================================

# Note that this function requires the packages "lme4" and "Matrix" to be installed


get_neff_mis <- function(model, N, t.points, survival) {
  
  n <- length(t.points)
  
  # extraction of model components
  reVar <- lme4::VarCorr(model)
  sigma2 <- sigma(model)^2
  Zt <- as.matrix(lme4::getME(model, "Zt")[1:2, 1:n, drop = FALSE])
  Z <- t(Zt)
  
  # pre-allocate and construct design matrices 
  ones <- rep(1, n)
  zeros <- rep(0, n)
  X_control <- cbind(intercept = ones, 
                     t = t.points, 
                     treat = zeros, 
                     `t:treat` = zeros)
  X_treat <- cbind(intercept = ones,
                   t = t.points,
                   treat = ones,
                   `t:treat` = t.points)
  
  # construct D, the variance covariance matrix of random effects
  D <- Matrix::bdiag(lapply(reVar, function(x) as(Matrix::bdiag(x), "generalMatrix")))
  
  # precompute all possible subsets under attrition
  row_indices <- lapply(n:1, function(k) if(k == n) 1:n else 1:k)
  
  # computation of V and W matrices
  V_subsets <- lapply(row_indices, function(rows) {
    Z_sub <- Z[rows, , drop = FALSE]
    as(Z_sub %*% D %*% t(Z_sub) + sigma2 * Matrix::Diagonal(length(rows)), "dgCMatrix")
  })
  
  W_subsets <- lapply(V_subsets, function(V) {
    Matrix::Diagonal(x = diag(as.matrix(V)))
  })
  
  # computation of inverses
  inverses <- lapply(seq_along(V_subsets), function(i) {
    list(
      V_inv = solve(V_subsets[[i]]),
      W_inv = solve(W_subsets[[i]])
    )
  })
  
  # empty objects for summed matrices
  sum_mat <- sum_mat_indep <- matrix(0, 4, 4)
  
  # backwards loop: from n (full data) to 1 (only first time point remains)
  # k represents the last time point of measurement
  for (k in n:1) {
    idx <- n - k + 1
    rows <- row_indices[[idx]]
    
    # use precomputed inverses
    inv <- inverses[[idx]]
    
    # accumulate sums
    X_control_sub <- X_control[rows, , drop = FALSE]
    X_treat_sub <- X_treat[rows, , drop = FALSE]
    
    sum_mat <- sum_mat + survival[k] * (
      t(X_control_sub) %*% inv$V_inv %*% X_control_sub +
        t(X_treat_sub) %*% inv$V_inv %*% X_treat_sub
    )
    
    sum_mat_indep <- sum_mat_indep + survival[k] * (
      t(X_control_sub) %*% inv$W_inv %*% X_control_sub +
        t(X_treat_sub) %*% inv$W_inv %*% X_treat_sub
    )
  }
  
  # matrix solving with fallback to pseudoinverse in case of singularity
  var_beta_hat <- tryCatch(
    solve(sum_mat),
    error = function(e) MASS::ginv(sum_mat)
  )
  
  var_betahat_indep <- tryCatch(
    solve(sum_mat_indep),
    error = function(e) MASS::ginv(sum_mat_indep)
  )
  
  # final computation of N_eff for the interaction coefficient
  w <- var_betahat_indep[4,4] / var_beta_hat[4,4]
  N_eff <- w * N * n
  
  return(N_eff)
}
