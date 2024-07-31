# METHOD COMPARISON



# Monte Carlo approach with increased samples
library(MASS)
library(mvtnorm)

# Monte Carlo approach with prior variance adjustment
calc_BF_monte_carlo <- function(N, mean_prior, hypothesis, mean_posterior, sigma_posterior, n_samples = 10000) {
  constraints <- unlist(strsplit(hypothesis, " "))
  params <- constraints[seq(1, length(constraints), 2)]
  
  mean_prior <- mean_prior[params]
  mean_posterior <- mean_posterior[params]
  
  region <- function(x) {
    names(x) <- names(mean_prior)
    result <- TRUE
    for (i in seq(1, length(constraints) - 2, 2)) {
      param1 <- constraints[i]
      constraint <- constraints[i + 1]
      param2 <- constraints[i + 2]
      
      if (constraint == ">") {
        result <- result && (x[param1] > x[param2])
      } else if (constraint == "<") {
        result <- result && (x[param1] < x[param2])
      } else if (constraint == "=") {
        result <- result && (x[param1] == x[param2])
      }
      if (!result) {
        return(FALSE)
      }
    }
    return(result)
  }
  
  samples_prior <- mvrnorm(n_samples, mu = mean_prior, Sigma = sigma_posterior * N)
  colnames(samples_prior) <- names(mean_prior)
  in_region_prior <- apply(samples_prior, 1, region)
  complexity <- mean(in_region_prior)
  
  samples_posterior <- mvrnorm(n_samples, mu = mean_posterior, Sigma = sigma_posterior)
  colnames(samples_posterior) <- names(mean_posterior)
  in_region_posterior <- apply(samples_posterior, 1, region)
  fit <- mean(in_region_posterior)
  
  BFu <- fit / complexity
  
  return(BFu)
}


# Numerical integration approach with prior variance adjustment
calc_BF_direct <- function(N, mean_prior, hypothesis, mean_posterior, sigma_posterior) {
  constraints <- unlist(strsplit(hypothesis, " "))
  params <- constraints[seq(1, length(constraints), 2)]
  
  mean_prior <- mean_prior[params]
  mean_posterior <- mean_posterior[params]
  
  region <- function(mean, sigma) {
    lower <- rep(-Inf, length(mean))
    upper <- rep(Inf, length(mean))
    
    for (i in seq(1, length(constraints) - 2, 2)) {
      param1 <- constraints[i]
      constraint <- constraints[i + 1]
      param2 <- constraints[i + 2]
      
      idx1 <- which(names(mean) == param1)
      idx2 <- which(names(mean) == param2)
      
      if (constraint == ">") {
        lower[idx1] <- mean[idx2]
      } else if (constraint == "<") {
        upper[idx1] <- mean[idx2]
      } else if (constraint == "=") {
        lower[idx1] <- mean[idx2]
        upper[idx1] <- mean[idx2]
      }
    }
    
    prob <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma, keepAttr = F)
    return(prob)
  }
  
  
  complexity <- region(mean=mean_prior, sigma=(sigma_posterior * N))
  fit <- region(mean=mean_posterior, sigma=sigma_posterior)
  
  BFu <- fit / complexity
  
  return(BFu)
}

# Numerical integration approach with prior variance adjustment by hand
calc_BF_direct_hand <- function(N, sigma, est, hypothesis, fraction=1) {
  b <- fraction/N
  struc <- unlist(strsplit(hypothesis, " "))
  params <- struc[seq(1, length(struc), 2)]
  constr <- struc[seq(2, length(struc)-1, 2)]

  # Contrast matrix
  Tmat <- matrix(0,2,3)
  
  if(all(constr==">")){
    lower <- rep(0, length(constr))
    upper <- rep(Inf, length(constr))
    Tmat[1,c(1,2)] <- Tmat[2,c(2,3)] <- c(1,-1)
  } else if(all(constr=="<")){
    lower <- rep(-Inf, length(constr))
    upper <- rep(0, length(constr))
    Tmat[1,c(1,2)] <- Tmat[2,c(2,3)] <- c(-1,1)
  } else {stop("not all constraints equal")}

  c <- mvtnorm::pmvnorm(lower=lower, upper=upper,
                        mean=rep(0, length(constr)), 
                        sigma=(Tmat %*% sigma %*% t(Tmat))/b,
                        keepAttr = F)
  f <- mvtnorm::pmvnorm(lower=lower,
                        upper=upper,
                        mean=c(Tmat %*% est),
                        sigma=Tmat %*% sigma %*% t(Tmat),
                        keepAttr = F) 
 BFu <- f/c
 return(BFu)
}


# Example usage
mean_prior <- c(a = 0, b = 0, c = 0)

mean_posterior <- est <- c(a = 3, b = 2, c = 1)
sigma_posterior <- matrix(c(1, 0.5, 0.5,
                            0.5, 1, 0.5,
                            0.5, 0.5, 1), nrow = 3)

sigma <- sigma_posterior

hypothesis <- "a > b > c"
N <- 100

BF_monte_carlo <- calc_BF_monte_carlo(N, mean_prior, sigma_prior, hypothesis, mean_posterior, sigma_posterior)
BF_direct <- calc_BF_direct(N, mean_prior, sigma_prior, hypothesis, mean_posterior, sigma_posterior)



BF_hand <- calc_BF_direct_hand(N=N, sigma=sigma, est=mean_posterior, hypothesis = hypothesis)
BF_pack <- BF(mean_posterior,
              Sigma=sigma,
              n=N,
              hypothesis ="a>b>c")
BF_pack <- BF_pack[["BFtable_confirmatory"]]


BF_monte_carlo
BF_direct
BF_hand
BF_pack

# SIMULATION 1 - growing effect
nsim <- 100

# Initialize a list to store the vectors
mean_post_list <- list()

# Loop to generate each vector
for (i in 0:(nsim-1)) {
  # Calculate the difference
  diff <- 2 * i / (nsim - 1)
  vec <- c(0, diff, 2 * diff)
  names(vec) <- c("a", "b", "c")
  mean_post_list[[i + 1]] <- vec
}

BFs_mc <- rep(NA, nsim)
BFs_dir <- rep(NA, nsim)
BFs_bain <- rep(NA, nsim)
BFs_pack <- rep(NA, nsim)
BFs_hand <- rep(NA, nsim)


for(i in 1:nsim){
  BFs_mc[i] <- calc_BF_monte_carlo(N=N, mean_prior=mean_prior,
                                   hypothesis=hypothesis,
                                   mean_posterior=mean_post_list[[i]],
                                   sigma_posterior=sigma_posterior)
  
  BFs_dir[i] <- calc_BF_direct(N=N, mean_prior=mean_prior,
                               hypothesis=hypothesis,
                               mean_posterior=mean_post_list[[i]],
                               sigma_posterior=sigma_posterior)
  
  a <- bain(n=N, hypothesis = hypothesis, Sigma=sigma_posterior, x=mean_post_list[[i]])
  BFs_bain[i] <- a[["fit"]][["BF.u"]][1]
  
  b <- BF(mean_post_list[[i]], Sigma=sigma, n=N, hypothesis = hypothesis)
  BFs_pack[i] <- b[["BFtable_confirmatory"]][1,6]
  
  BFs_hand[i] <- calc_BF_direct_hand(N=N, sigma=sigma, est=mean_post_list[[i]], hypothesis=hypothesis)
  
  print(i/nsim)
}

plot(x=seq(1:nsim), y=BFs_hand, type="l")
lines(x=BFs_dir, col="blue")
lines(BFs_bain, col="red")
lines(x=BFs_mc, col="orange")
lines(x=BFs_pack, col="green")


# SIMULATION 2 - only one difference becomes larger

nsim2 <- 100
mean_prior <- c(a = 0, b = 0, c = 0)
sigma_prior <- matrix(c(1, 0.3, 0.3,
                        0.3, 1, 0.3,
                        0.3, 0.3, 1), nrow = 3)

mean_posterior <- c(a = 0, b = 0.2, c = 0.4)
sigma_posterior <- matrix(c(1, 0.5, 0.5,
                            0.5, 1, 0.5,
                            0.5, 0.5, 1), nrow = 3)

N <- seq(10, 10000, length.out=nsim2)
hypothesis <- "a < b < c"

BFs_mc2 <- rep(NA, nsim)
BFs_dir2 <- rep(NA, nsim)
BFS_bain2 <- rep(NA, nsim)


for(i in 1:nsim2){
  BFs_mc2[i] <- calc_BF_monte_carlo(N=N[i], mean_prior=mean_prior,
                                   sigma_prior=sigma_prior,
                                   hypothesis=hypothesis,
                                   mean_posterior=mean_posterior,
                                   sigma_posterior=sigma_posterior)
  
  BFs_dir2[i] <- calc_BF_direct(N=N[i], mean_prior=mean_prior,
                               sigma_prior=sigma_prior,
                               hypothesis=hypothesis,
                               mean_posterior=mean_posterior,
                               sigma_posterior=sigma_posterior)
  
  a <- bain(x=mean_posterior, Sigma=sigma_posterior, hypothesis=hypothesis, n=N[i])
  BFs_bain2[i] <- a[["fit"]][["BF.u"]][1]
  print(i/nsim2)
}

plot(x=seq(1:nsim2), y=BFs_mc2, type="l")
lines(x=BFs_bain2, col="red")


