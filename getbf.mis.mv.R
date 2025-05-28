#===============================================================================
# Function for data generation and Bayes Factor calculation for a single data set
# with Missing Data Patterns for multiple parameters
#===============================================================================

# Hypotheses to be tested:
# H0: beta2 = 0
# H1: beta2 > 0
# where beta2 is the coefficient of interaction of time and treatment condition (treatment effect)

# The function "getbf" uses the following arguments:
# N = the total sample size (number of subjects)
# t.points = position of the measurement occasions in time
# var.u0 = intercept variance
# var.u1 = slope variance
# cov = covariance between intercept and slope variance
# var.e = error variance
# eff.size = effect size defined as beta/sqrt(var.u1), where beta is the coefficient of interaction
# fraction = fraction of information used to specify prior, b = fraction/N
# where n = number of measurement occasions
# log.grow = indicates whether to use logarithmic (TRUE) or linear growth (FALSE)
# hyp = for which hypotheses should the BF be calculated ("H0", "H1", "both" or "h0", "h1", "b")

# Note 1: this function requires the function "get_neff_mis" to be loaded into the global environment
# This function requires the packages "MASS", "dplyr", and "lme4" to be installed
#-------------------------------------------------------------------------------

library(bain)
library(dplyr)

getbf_mis_mv <- function(N, attrition="weibull", params=list(.5,1), hypothesis,
                         t.points, var.u0, var.u1, cov, var.e, eff.sizes, fraction, 
                         log.grow, beta1){
  
  if(attrition=="exponential" && length(params)!=1){stop("The exponential distribution requires only one parameter in the 'params' argument")}
  if((attrition=="weibull" | attrition=="linear-exponential" | attrition=="log-logistic" | attrition=="gompertz") && length(params)!=2){stop("The weibull, linear-exponential, log-logistic, and gompertz distribution require two parameters in the 'params' argument")}
  if(attrition=="nonparametric" && length(params)!=length(t.points)){stop("In case of a nonparametric survival function, the number of 'params' must equal the number of timepoints")}
  if(attrition=="nonparametric" && all(params == cummin(params))==F){stop("In case of a nonparametric survival function, the elements of 'params' must be monotonically decreasing")}
  
  n <- length(t.points)  # number of measurement occasions
  # create time variable t
  if(log.grow==F) {t <- rep(t.points, N) # if logarithmic growth is used, take log of t.points
  } else {if(min(t.points)==0) {t <- rep(log(t.points+1), N)} else {t <- rep(log(t.points), N)} # if the first timepoint is zero, we add 1 to all timepoints because log(0) is undefined
  } 
  t.prop <- t.points/max(t.points) # create rescaled time variable (from 0 to 1)
  id <- rep(seq_len(N), each=n)  # create ID variable
  treat <- as.character(gl(n=3, k=n, length=N*n, labels=c("a","b","c"))) # create treatment variable
  dat0 <- data.frame(id, treat, t) # combine into data frame
  dat0$treat <- factor(dat0$treat, levels = c("a", "b", "c")) # Forces "a" as reference
  multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, cov, cov, var.u1), nrow=2, ncol=2)) # draw random effects
  
  
  
  compute_survival <- function(distributions, params, t.points) {
    # Convert single distribution/params to list format
    if (is.character(distributions) && length(distributions) == 1) {
      distributions <- list(distributions)
    }
    if (is.numeric(params)) {
      params <- list(params)
    }
    
    # Validate inputs
    if (!is.list(params) || !all(sapply(params, is.numeric))) {
      stop("params must be a list of numeric vectors")
    }
    if (length(distributions) != length(params)) {
      stop("Number of distributions must match number of parameter sets")
    }
    
    # Convert time points to proportional time (0-1)
    t.prop <- t.points/max(t.points)
    
    # Define distribution functions (same as before)
    dist_functions <- list(
      weibull = function(pars, t) {
        if (length(pars) < 2) stop("Weibull needs omega, gamma")
        (1-pars[1])^(t^pars[2])
      },
      modified_weibull = function(pars, t) {
        if (length(pars) < 3) stop("Modified Weibull needs omega, gamma, kappa")
        (exp(time^pars[2]*exp(pars[3]*(time-1))*log(1-pars[1])))
      },
      exponential = function(pars, t) {
        if (length(pars) < 1) stop("Exponential needs omega")
        (1-pars[1])^t
      },
      log.logistic = function(pars, t) {
        if (length(pars) < 2) stop("Log-logistic needs omega, gamma")
        (1-pars[1])/((1-pars[1]) + pars[1]*t^pars[2])
      },
      gompertz = function(pars, t) {
        if (length(pars) < 2) stop("Gompertz needs omega, gamma")
        exp((log(1-pars[1])/(exp(pars[2])-1))*(exp(pars[2]*t)-1))
      },
      nonparametric = function(pars, t) {
        if (t[1] == 0) {
          if (length(pars) < max(t)+1) stop("Need ", max(t)+1, " params for t=0 start")
          pars
        } else {
          if (length(pars) < max(t)) stop("Need ", max(t), " params for t=1 start")
          pars
        }
      }
    )
    
    # Calculate all survival curves
    results <- Map(function(dist, pars) {
      if (!dist %in% names(dist_functions)) {
        stop("Unknown distribution: ", dist)
      }
      dist_functions[[dist]](pars, t.prop)
    }, distributions, params)
    
    # Return vector if single distribution, else list
    if (length(results) == 1) {unlist(results)} else {results}
  }
  
  # Define distributions and parameters for 2 groups
  attrition <- c("weibull", "gompertz")
  params <- list(c(0.5, 1), c(.5, .5))

  # Compute survival curves
  if(any(attrition != F)){
    survival_curves <- compute_survival(attrition, params, t.points)
    # hazard_curves <- lapply(survival_curves, function())
  }


  # create missingness with survival function 
  # if(dropout == TRUE){
  #   if (distribution == "weibull"){                                    # Weibull
  #     weibull <- function(omega, gamma, time){ 
  #       (1-omega)^time^gamma
  #     }
  #     survival <- weibull(omega=params[[1]], gamma=params[[2]], time=t.prop)
  #   } else if (distribution == "exponential"){                     # Exponential
  #     exponential <- function(omega, time){
  #       (1-omega)^time
  #     }
  #     survival <- exponential(omega=params[[1]], time=t.prop)
  #   } else if (distribution == "log-logistic"){                   # Log-logistic
  #     log_logistic <- function(omega, gamma, time){
  #       (1-omega)/((1-omega)+omega*time^gamma)
  #     }
  #     survival <- log_logistic(omega=params[[1]], gamma=params[[2]], time=t.prop)
  #   } else if (distribution == "linear-exponential"){       # Linear-exponential
  #     linear_exponential <- function(omega, gamma, time){
  #       (exp((.5*gamma+log(1-omega))*time-.5*gamma*time^2))
  #     }
  #     survival <- linear_exponential(omega=params[[1]], gamma=params[[2]], time=t.prop)
  #   } else if(distribution == "modified_weibull"){            # Modified Weibull
  #     mod_weibull <- function(omega, gamma, kappa, time){
  #       (exp(time^gamma*exp(kappa*(time-1))*log(1-omega)))
  #     }
  #     survival <- mod_weibull(omega=params[[1]], gamma=params[[2]], kappa=params[[3]], time=t.prop)
  #   } else if(distribution=="gompertz"){                              # Gompertz
  #     gompertz <- function(omega, gamma, time){
  #       (exp((log(1-omega)/(exp(gamma)-1))*(exp(gamma*time)-1)))
  #     }
  #     survival <- gompertz(omega=params[[1]], gamma=params[[2]], time=t.prop)
  #   } else if(distribution=="nonparametric"){                   # non-parametric
  #     survival <- unlist(params)
  #   }
  #   
    # compute hazard: (S_t - S_{t+1}) / S_t
    shifted_survival <- c(survival[-1], NA)  # drop first element and append NA
    hazard <- (survival - shifted_survival) / survival  
  } else if(attrition==F){
    survival <- rep(1, n)
  }

  # generate data under the research hypothesis
  u0 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n) # random intercepts
  u1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n) # random slopes
  e <- rnorm(N*n, 0, sqrt(var.e)) # error variance for H1
  
  # Create treatment contrast variables for modeling
  treat_TAU <- as.numeric(treat == "b")
  treat_INT <- as.numeric(treat == "c")
  
  beta2_TAU <- unlist(eff.sizes[1]) * sqrt(var.u1) # create coefficient of interaction for b from effect size
  beta3_INT <- unlist(eff.sizes[2]) * sqrt(var.u1) + beta2_TAU # create coefficient of interaction for c from effect size
  
  y <- beta1*t + u0 + beta2_TAU*treat_TAU*t + beta3_INT*treat_INT*t + u1*t + e # create outcome variable y under H1
  if(attrition == T){
    dat <- data.frame(dat0, y, hazard)
    dat[dat$t.prop==1, "hazard"] <- NA # h(t_last) is undefined
    suppressWarnings(dat$mis <- rbinom(n=nrow(dat), size=1, prob=dat$hazard)) # suppress warning about NAs being produced
    dat <- dat %>% group_by(id) %>% mutate(mis = ifelse(cumany(mis == 1), 1, mis))
    dat$y[which(dat$mis==1)] <- NA
  } else if (attrition == F){
    dat <- data.frame(dat0, y)}

  # fit MLM to dataset
  models <- lme4::lmer(y ~ t + treat + t:treat + (1 + t | id), data = dat, REML=F, control = lme4::lmerControl(calc.derivs = F))  # fit MLM model under H1
  est <- models@beta[c(2,5,6)] # extract estimates of beta2 and beta3 under H0
  names(est) <- c("a", "b", "c")
  sig_WL <- as.matrix(vcov(models)[2,2])  # extract variance of estimates under H0  
  sig_TAU <- as.matrix(vcov(models)[5,5])  # extract variance of estimates under H0  
  sig_INT <- as.matrix(vcov(models)[6,6])  # extract variance of estimates under H0  
  
  # calculate N_eff
  n_eff <- get_neff_mis_mv(model=models, N=N, t.points=t.points, survival=survival)

  bf_res <- bain(x=est, Sigma=list(sig_WL, sig_TAU, sig_INT), n=n_eff, 
                 hypothesis=hypothesis, group_parameters = 1, joint_parameters = 0)
  
  # extract number of different hypotheses to be evaluated
  n_hyp <- length(gregexpr(";", hypothesis, fixed = TRUE)[[1]]) + 1
  
  bf_c <- bf_res[["fit"]][["BF.c"]][1]
  bfs <- bf_res[["BFmatrix"]]
  PMPc <- bf_res[["fit"]][["PMPc"]][1]

  return(output = list(bf_c=bf_c, PMPc=PMPc))
}
