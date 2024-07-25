# Function to calculate the BFu analytically for informative hypotheses
# bf.hand requires the function "hyp2mat" to be loaded into environment 
# WARNING: bf.hand cannot yet handle hypotheses with heterogeneous constraints!

bf.hand <- function(N, est, sigma, hypothesis, fraction=1) {
  
  # Add spaces around inequality signs if they are not present
  hypothesis <- gsub("([>=<])", " \\1 ", hypothesis)
  # Remove any extra spaces that may have been introduced
  hypothesis <- gsub("\\s+", " ", hypothesis)
  # Trim leading and trailing spaces
  hypothesis <- trimws(hypothesis)
  
  b <- fraction/N
  struc <- unlist(strsplit(hypothesis, " "))
  params <- struc[seq(1, length(struc), 2)]
  constr <- struc[seq(2, length(struc)-1, 2)]
  
  # Contrast matrix
  Rmat <- hyp2mat(hypothesis = hypothesis)
  Rmat.ineq <- Rmat$Ineq.mat
  Rmat.eq <- Rmat$Eq.mat
  
  # Define bounds for inequalities
  lower <- rep(-Inf, length(constr))
  upper <- rep(Inf, length(constr))
  
  # Assign bounds based on inequality constraints
  if(any(constr == ">")) {
    lower[constr == ">"] <- 0
  }
  if(any(constr == "<")) {
    upper[constr == "<"] <- 0
  }
  
  # Compute the complexity (c) for inequalities
  c <- 1
  if(!is.null(Rmat.ineq) && nrow(Rmat.ineq) > 0) {
    c <- mvtnorm::pmvnorm(lower=lower, upper=upper,
                          mean=rep(0, length(constr)), 
                          sigma=(Rmat.ineq %*% sigma %*% t(Rmat.ineq))/b,
                          keepAttr = F)
  }
  
  # Compute the fit (f) for inequalities
  f <- 1
  if(!is.null(Rmat.ineq) && nrow(Rmat.ineq) > 0) {
    if(any(constr == "<")) {
      est_adjusted <- -est
    } else {
      est_adjusted <- est
    }
    f <- mvtnorm::pmvnorm(lower=lower,
                          upper=upper,
                          mean=c(Rmat.ineq %*% est_adjusted),
                          sigma=Rmat.ineq %*% sigma %*% t(Rmat.ineq),
                          keepAttr = F) 
  }
  
  # Handle the equality constraints
  if(!is.null(Rmat.eq) && nrow(Rmat.eq) > 0) {
    sigma_eq <- Rmat.eq %*% sigma %*% t(Rmat.eq)
    f <- mvtnorm::dmvnorm(x=rep(0, length(constr)), mean=Rmat.eq %*% est, sigma=sigma_eq, log=F)
    c <- mvtnorm::dmvnorm(rep(0, nrow(Rmat.eq)), sigma=sigma_eq/b, log=F)
  }
  
  BFu <- f/c
  return(list(complexity=c, fit=f, BFu=BFu))
}


# example values
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

# results are the same, yay!
