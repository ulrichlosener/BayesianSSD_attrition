
bf_hand <- function(N, est, sigma, hypothesis, fraction = 1) {
  
  hypothesis <- gsub("([>=<])", " \\1 ", hypothesis) # Add spaces
  hypothesis <- gsub("\\s+", " ", hypothesis) # Remove any extra spaces
  hypothesis <- trimws(hypothesis) # Trim leading and trailing spaces

  b <- fraction/N
  struc <- unlist(strsplit(hypothesis, " "))
  params <- struc[seq(1, length(struc), 2)]
  constr <- struc[seq(2, length(struc) - 1, 2)]
  ineq.constr <- constr[constr == "<" | constr == ">"]
  eq.constr <- constr[constr == "="]
  
  # Contrast matrix
  Rmat <- hyp2mat(hypothesis = hypothesis)
  Rmat_ineq <- Rmat$Ineq.mat
  Rmat_eq <- Rmat$Eq.mat
  
  # Compute the complexity (c) and fit (f) for inequalities
  if(!is.null(Rmat_ineq) && nrow(Rmat_ineq) > 0) {
    
    if(any(constr == "<")) {
      est_adjusted <- -est
    } else {
      est_adjusted <- est
    }
    
    # Initialize bounds for inequalities
    lower <- rep(-Inf, length(ineq.constr))
    upper <- rep(Inf, length(ineq.constr))
    
    # Assign bounds based on inequality constraints
    lower[ineq.constr == ">"] <- 0
    upper[ineq.constr == "<"] <- 0
    
    c_ineq <- as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                                          mean = rep(0, length(ineq.constr)), 
                                          sigma = Rmat_ineq %*% sigma %*% t(Rmat_ineq)/b,
                                          keepAttr = F))
    f_ineq <- as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                                          mean = c(Rmat_ineq %*% est_adjusted),
                                          sigma = Rmat_ineq %*% sigma %*% t(Rmat_ineq),
                                          keepAttr = F))
    BFu_ineq <- f_ineq / c_ineq
  } else {
    c_ineq <- 1
    f_ineq <- 1
    BFu_ineq <- 1
  }
  
  if (!is.null(Rmat_eq) && nrow(Rmat_eq) > 0) {
    c_eq <- mvtnorm::dmvnorm(rep(0, nrow(Rmat_eq)), 
                             sigma = Rmat_eq %*% sigma %*% t(Rmat_eq)/b)
    f_eq <- mvtnorm::dmvnorm(x = rep(0, nrow(Rmat_eq)), mean = Rmat_eq %*% est, 
                             sigma = Rmat_eq %*% sigma %*% t(Rmat_eq))
    BFu_eq <- f_eq / c_eq
  } else {
    c_eq <- 1
    f_eq <- 1
    BFu_eq <- 1
  }
  
  BFu <- BFu_ineq * BFu_eq
  
  return(list(complex_ineq = c_ineq, 
              complex_eq = c_eq,
              fit_ineq = f_ineq, 
              fit_eq = f_eq,
              BFu = BFu))
}

N <- 100
est <- c(1,2,3)
names(est) <- c("a", "b", "c")
sigma <- matrix(c(1, 0.3, 0.3,
                  0.3, 1, 0.3,
                  0.3, 0.3, 1), nrow = 3)
hypothesis <- "a=b<c"

bf_hand2(N=100, est=est, sigma=sigma, hypothesis=hypothesis)



a <- BF(est, hypothesis=hypothesis, Sigma = sigma, n=N)
a[["BFtable_confirmatory"]]

bain(x=est, hypothesis=hypothesis, Sigma=sigma, n=N)





# emore xample values
N <- 100
est <- c(3,2,1)
names(est) <- c("a", "b", "c")
est3 <- c(3,2,1,0)
names(est3) <- c("a", "b", "c", "d")
sigma <- matrix(c(1, 0.3, 0.3,
                  0.3, 1, 0.3,
                  0.3, 0.3, 1), nrow = 3)
sigma3 <- matrix(c(1, 0.3, 0.2, 0.1,
                   0.3, 1, 0.3, 0.2,
                   0.2, 0.3, 1, 0.3,
                   0.1, 0.2, 0.3, 1), nrow = 4)
hypothesis1 <- "a= b = c"
hypothesis2a <- "a> b > c"
hypothesis2b <- "a < b< c"
hypothesis3 <- "a > b>c>d"

# compare Bf.hand with BF from BFpack for hypothesis 1
bf.hand(N=N, sigma=sigma, est=est, hypothesis=hypothesis1)
a1 <- BF(x=est, Sigma = sigma, n=N, hypothesis = hypothesis1)
a1[["BFtable_confirmatory"]][1,c(1,3,5)]

# compare Bf.hand with BF from BFpack for hypothesis 2a
bf.hand(N=N, sigma=sigma, est=est, hypothesis=hypothesis2a)
a4 <- BF(x=est, Sigma = sigma, n=N, hypothesis = hypothesis2a)
a4[["BFtable_confirmatory"]][1,c(2,4,6)]

# compare Bf.hand with BF from BFpack for hypothesis 2b
bf.hand(N=N, sigma=sigma, est=est, hypothesis=hypothesis2b)
a5 <- BF(x=est, Sigma = sigma, n=N, hypothesis = hypothesis2b)
a5[["BFtable_confirmatory"]][1,c(2,4,6)]


# compare Bf.hand with BF from BFpack for hypothesis 3
bf.hand(N=N, sigma=sigma3, est=est3, hypothesis=hypothesis3)
a3 <- BF(x=est3, Sigma = sigma3, n=N, hypothesis = hypothesis3)
a3[["BFtable_confirmatory"]][1,c(2,4,6)]
