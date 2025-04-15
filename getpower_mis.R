#===============================================================================
# Function for power calculation with Missing Data Patterns
#===============================================================================


getpower_mis <- function(dropout=T, omega=.3, gamma=1, m=1000, N=72, t.points=c(0,1,2,3,4), 
                         var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
                         BFthres=3, fraction=1, beta1=0, log.grow=F, seed=NULL, hyp="both",
                         test="alt"){
  
  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
  
  plan(multisession, workers = availableCores() - 1)  # Use all but one core
  
  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  suppressMessages({
    bfs <- future_lapply(Ns, getbf_mis, dropout=dropout, omega=omega, gamma=gamma, 
                         t.points=t.points, var.u0=var.u0, var.u1=var.u1, 
                         cov=cov, var.e=var.e, eff.size=eff.size, fraction=fraction, 
                         log.grow=log.grow, hyp=hyp, beta1=beta1, future.seed=TRUE)
  })
  
  plan(sequential)  # Reset plan to avoid unexpected parallel behavior later
  
  
  # comparison of H0 and H1 ----------------------------------------------------
  if((hyp == "both" | hyp == "b") & test == "alt"){
    bf0 <- sapply(bfs, "[[", 1) # extract all BF01
    bf1 <- sapply(bfs, "[[", 2) # extract all BF10
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0[bf0>=BFthres])/m,
                power.H1 = length(bf1[bf1>=BFthres])/m))
    
  } else if((hyp == "h0" | hyp == "H0") & test == "alt"){
    bf0 <- sapply(bfs, "[[", 1)
    
    # return power for H0 vs H1 only
    return(list(power.H0 = length(bf0[bf0>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "1") & test == "alt"){
    bf1 <- sapply(bfs, "[[", 1)
    
    # return power for H1 vs H0 only
    return(list(power.H1 = length(bf1[bf1>=BFthres])/m))
    
    # comparison of H and Hu ---------------------------------------------------
  } else if((hyp == "both" | hyp == "b") & (test == "hu" | test == "Hu")){
    bf0u <- sapply(bfs, "[[", 3) # extract all BF0u
    bf1u <- sapply(bfs, "[[", 4) # extract all BF1u
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0u[bf0u>=BFthres])/m,
                power.H1 = length(bf1u[bf1u>=BFthres])/m)) 
    
  } else if((hyp == "h0" | hyp == "H0") & (test == "hu" | test == "Hu")){
    bf0u <- sapply(bfs, "[[", 2)
    
    # return power for H0 vs Hu only
    return(list(power.H0 = length(bf0u[bf0u>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "H1") & (test == "hu" | test == "Hu")){
    bf1u <- sapply(bfs, "[[", 2)
    
    # return power for H1 vs Hu only
    return(list(power.H1 = length(bf1u[bf1u>=BFthres])/m))
    
    # comparison of H and Hc ---------------------------------------------------
  } else if((hyp == "both" | hyp == "b") & (test == "hc" | test == "Hc")){
    bf0c <- sapply(bfs, "[[", 5) # extract all BF0u
    bf1c <- sapply(bfs, "[[", 6) # extract all BF1u
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0c[bf0c>=BFthres])/m,
                power.H1 = length(bf1c[bf1c>=BFthres])/m)) 
    
  } else if((hyp == "h0" | hyp == "H0") & (test == "hc" | test == "Hc")){
    bf0c <- sapply(bfs, "[[", 3)
    
    # return power for H0 vs Hc only
    return(list(power.H0 = length(bf0c[bf0c>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "H1") & (test == "hc" | test == "Hc")){
    bf1c <- sapply(bfs, "[[", 3)
    
    # return power for H1 vs Hc only
    return(list(power.H1 = length(bf1c[bf1c>=BFthres])/m))
  }
}

# END OF FUNCTION --------------------------------------------------------------
