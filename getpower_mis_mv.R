#===============================================================================
# Function for power calculation with Missing Data Patterns for three groups
#===============================================================================

library(future)
library(future.apply)

getpower_mis_mv <- function(attrition="weibull", params=c(.5,1), 
                         m=10, N=100, t.points=c(0,1,2,3,4), var.u0=0.03, 
                         var.u1=.1, var.e=.02, cov=0, eff.sizes=c(0, .5, .5), 
                         BFthres=5, fraction=1, log.grow=F, seed=NULL, 
                         hypothesis="a<b<c", PMPthres=.9){
  
  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
  
  plan(multisession, workers = availableCores() - 1)  # Use all but one core
  
  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  suppressMessages({
    bfs <- future_lapply(Ns, getbf_mis_mv, attrition=attrition, 
                         params=params, hypothesis=hypothesis, t.points=t.points, 
                         var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, 
                         eff.sizes=eff.sizes, fraction=fraction, log.grow=log.grow, 
                         future.seed = TRUE)
  })
  
  plan(sequential)  # Reset plan to avoid unexpected parallel behavior later
  
  bf <- sapply(bfs, function(x) x[1])
  pmp <- sapply(bfs, function(x) x[2])

  power_bf <- mean(bf > BFthres)
  power_pmp <- mean(pmp > PMPthres)
  
  return(list(power_bf=power_bf,
              power_pmp=power_pmp))
}

# END OF FUNCTION --------------------------------------------------------------
