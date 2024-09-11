bf_hand2 <- function(N, est, sigma, hypothesis, fraction = 1) {
  
  # Clean and prepare the hypothesis string
  hypothesis <- gsub("([>=<])", " \\1 ", hypothesis) # Add spaces around inequality operators
  hypothesis <- gsub("\\s+", " ", hypothesis) # Remove any extra spaces
  hypothesis <- trimws(hypothesis) # Trim leading and trailing spaces
  
  b <- fraction/N
  struc <- unlist(strsplit(hypothesis, " "))
  params <- struc[seq(1, length(struc), 2)]
  constr <- struc[seq(2, length(struc) - 1, 2)]
  ineq.constr <- constr[constr == "<" | constr == ">"]
  eq.constr <- constr[constr == "="]
  
  # Generate the contrast matrix based on the hypothesis
  Rmat <- hyp2mat(hypothesis = hypothesis)
  Rmat_ineq <- Rmat$Ineq.mat
  Rmat_eq <- Rmat$Eq.mat
  
  # Compute the complexity (c) and fit (f) for inequalities
  if (!is.null(Rmat_ineq) && nrow(Rmat_ineq) > 0) {
    
    # Initialize bounds for inequalities (always [0, Inf])
    lower <- rep(0, length(ineq.constr))
    upper <- rep(Inf, length(ineq.constr))
    
    # Transform estimates and adjust sign for "<" constraints
    transformed_est <- c(Rmat_ineq %*% est)  # Transformation of the estimate
    
    # Compute complexity for inequality constraints
    c_ineq <- as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                                          mean = rep(0, length(ineq.constr)), 
                                          sigma = Rmat_ineq %*% sigma %*% t(Rmat_ineq) / b,
                                          keepAttr = F))
    
    # Compute fit for inequality constraints
    f_ineq <- as.numeric(mvtnorm::pmvnorm(lower = rep(0, length(ineq.constr)), upper = rep(Inf, length(ineq.constr)),
                                          mean = transformed_est,  # Adjusted estimate
                                          sigma = Rmat_ineq %*% sigma %*% t(Rmat_ineq),
                                          keepAttr = F))
    
    BFu_ineq <- f_ineq / c_ineq
  } else {
    c_ineq <- 1
    f_ineq <- 1
    BFu_ineq <- 1
  }
  
  # Compute the complexity (c) and fit (f) for equalities
  if (!is.null(Rmat_eq) && nrow(Rmat_eq) > 0) {
    c_eq <- mvtnorm::dmvnorm(rep(0, nrow(Rmat_eq)), 
                             sigma = Rmat_eq %*% sigma %*% t(Rmat_eq) / b)
    f_eq <- mvtnorm::dmvnorm(x = rep(0, nrow(Rmat_eq)), mean = Rmat_eq %*% est, 
                             sigma = Rmat_eq %*% sigma %*% t(Rmat_eq))
    BFu_eq <- f_eq / c_eq
  } else {
    c_eq <- 1
    f_eq <- 1
    BFu_eq <- 1
  }
  
  # Compute the final Bayes factor
  BFu <- BFu_ineq * BFu_eq
  
  return(list(complex_ineq = c_ineq, 
              complex_eq = c_eq,
              fit_ineq = f_ineq, 
              fit_eq = f_eq,
              BFu = BFu))
}
