#===============================================================================
# Function to calculate survival curves to model attrition
#===============================================================================

# "distributions" can be either a vector of length 1 or of length g, where g is the number of treatment conditions. 
# "params" is a list with either 1 elemtent or g elements. Each element contains the parameters for the respective survival function.
# "t.points" is a vector containing the measurement occasions.


survival <- function(distributions, params, t.points) {
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
  
  # Convert time points to proportional time (0-1)
  time <- t.points/max(t.points)
  
  # Define distribution functions
  dist_functions <- list(
    weibull = function(pars, t) {
      if (length(pars) < 2) stop("Weibull needs omega, gamma")
      (1-pars[1])^(t^pars[2])
    },
    modified_weibull = function(pars, t) {
      if (length(pars) < 3) stop("Modified Weibull needs omega, gamma, kappa")
      (exp(time^pars[2]*exp(pars[3]*(time-1))*log(1-pars[1])))
    },
    linear_exponential = function(pars, t) {
      if (length(pars) < 2) stop("Linear-exponential needs omega, gamma")
      # (exp((.5*pars[2]+log(1-pars[1]))*time - .5*pars[2]*time^2))
      (exp((.5*pars[2] + log(1-pars[1]))*time - (.5*pars[2]*time)))
    },
    loglogistic = function(pars, t) {
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
  
  # Calculate survival curve(s)
  results <- Map(function(dist, pars) {
    if (!dist %in% names(dist_functions)) {
      stop("Unknown distribution: ", dist)
    }
    dist_functions[[dist]](pars, time)
  }, distributions, params)
  
  return(results)
}
