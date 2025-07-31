
BayeSSD <- function(eta=.8, attrition="weibull", params=c(.5,1), 
                        m=100, t.points=c(0,1,2,3,4), var.u0=0.01, 
                        var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .8, .8), 
                        BFthres=5, fraction=1, log.grow=F, seed=NULL, 
                        hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
                        N_max=1000) {
  
  # error and warning messages in case of incorrect input
  if(eta<0 | eta>1) {stop("'eta' (the desired power level) must be between 0 and 1.")}
  if(m%%1!=0 | m<1) {stop("'m' must be a positive integer.")}
  if(!is.logical(log.grow)) {stop("'log.grow' must be either TRUE or FALSE.")}
  if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE.")}
  if(any(t.points<0)) {stop("all time points must be positive.")}
  if(var.u0<0 | var.u1<0 | cov<0 | var.e<0) {stop("all variance components must be positive.")}
  if(BFthres<0) {stop("'BFthres' must be positive.")}
  if(fraction%%1!=0 | fraction<1) {stop("'fraction' must be a positive integer, b=fraction/N.")}
  if(m<1000) {warning("Results with less than 1000 generated datasets per iteration can be unreliable.")}

  start_time <- Sys.time()

  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
 
  N <- list()
  
  Nmin <- 30            # (initial) minimal sample size 
  Nmax <- N_max         # (initial) maximum sample size
  condition <- FALSE    # condition initially FALSE until power criterion is reached
  j <- 1                # iteration counter
  
      while(condition == F){
        
        N[j] <- round((Nmin + Nmax)/2 - .1, digits = 0)  # current N is the mid point between Nmin and Nmax, rounded to the lower number
        # generate data and store BFs
        results <- getpower_mis_mv(attrition=attrition, params=params, m=m, N=unlist(N[j]), 
                                   log.grow=log.grow, fraction=fraction, 
                                    t.points=t.points, var.u0=var.u0, var.u1=var.u1, 
                                    cov=cov, var.e=var.e, eff.sizes=eff.sizes, 
                                    BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

        # check if condition is met
        ifelse(results$power_bf>=eta, 
               Nmax <- unlist(N[j]) - 1,
               Nmin <- unlist(N[j]) + 1
        )
        
        # Calculate time metrics
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        avg_time_per_iter <- elapsed / j
        remaining_time <- avg_time_per_iter * (16 - j)  # max_iter = 16
        
        # Print progress
        cat(
          sprintf("Iter %d: N = %d, Power = %.3f | Elapsed: %.1f minutes | Remaining: ~ %.1f minutes \n",
                  j, unlist(N[[j]]), results$power_bf, elapsed, remaining_time)
        )
        
        # if N increases by only 1 or f power level is very close to desired power level, condition is met and the algorithm stops
        if ((N[j] == Nmin+1 | Nmax == Nmin) | round(abs(results$power_bf - eta), 8) <= tol) {
          condition <- TRUE
          total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
          cat(sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f)\n",
                      j, total_time, unlist(N[[j]]), results$power_bf))
        }

        # increase iteration by 1
        j <- j+1
      }
  
  
}


# END OF FUNCTION --------------------------------------------------------------
