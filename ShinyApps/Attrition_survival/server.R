# Server part of the shiny app

server <- function(input, output) {

  
### Weibull function
 output$weibullplots <- renderPlot({
    
    validate(need(input$d_weib>=3,"Duration of the study should be at least three."))
    if(input$omega_weib<0 || input$omega_weib>1) {
      validate("`Omega` should be between zero and one.")
    }
    validate(need(input$gamma_weib>=0,"`Gamma` should be positive."))
    
    omega <- input$omega_weib
    gamma <- input$gamma_weib
    f <- input$f_weib
    D <- input$d_weib
    t.points <- seq(0, D, by=1/f)
    time <- seq(0, 1, length=length(t.points))
    remain <- (1-omega)^time^gamma
    n <- length(remain)
    
    hazard <- rep(NA, n)
    for(i in 1:(n-1)) {
      hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
    }
    
  par(mfrow=c(2,1))
  plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="probability", main="Survival Function")
  plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="Measurement Occasion", ylab="Probability", main="Hazard Function")
    }, height=600)
 
### Modified Weibull function
 output$mod_weibullplots <- renderPlot({
    
   validate(need(input$d_weib>=3,"Duration of the study should be at least three."))
   if(input$omega_mod_weib<0 || input$omega_mod_weib>1) {
     validate("`Omega` should be between zero and one.")
   }
   validate(need(input$gamma_mod_weib>=0, "`Gamma` should be positive."))
   validate(need(input$kappa_mod_weib>=0, "`Kappa` should be positive."))
   
   omega <- input$omega_mod_weib
   gamma <- input$gamma_mod_weib
   kappa <- input$kappa_mod_weib
   f <- input$f_mod_weib
   D <- input$d_mod_weib
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- exp(time^gamma*exp(kappa*(time-1))*log(1-omega))

   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="probability", main="Survival Function")
   plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="Measurement Occasion", ylab="Probability", main="Hazard Function")
 }, height=600)
 
 ### Linear-Exponential function
 output$lin_expplots <- renderPlot({
   
   validate(need(input$d_lin_exp>=3,"Duration of the study should be at least three."))
   if(input$omega_lin_exp < 0 || input$omega_lin_exp > 1) {
     validate("`Omega` should be between zero and one.")
   }
   if(input$gamma_lin_exp < -1 || input$gamma_lin_exp > 1) {
     validate("`Gamma` should be between minus one and one.")
   }

   omega <- input$omega_lin_exp
   gamma <- input$gamma_lin_exp
   kappa <- input$kappa_lin_exp
   f <- input$f_lin_exp
   D <- input$d_lin_exp
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- exp((.5*gamma + log(1-omega))*time - (.5*gamma*time^2))

   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="probability", main="Survival Function")
   plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="Measurement Occasion", ylab="Probability", main="Hazard Function")
 }, height=600)
 
 ### Log-Logistic function
 output$logplots <- renderPlot({
   
   validate(need(input$d_log>=3,"Duration of the study should be at least three."))
   if(input$omega_log < 0 || input$omega_log > 1) {
     validate("`Omega` should be between zero and one.")
   }
   validate(need(input$gamma_log>=0, "`Gamma` should be positive."))
   
   omega <- input$omega_log
   gamma <- input$gamma_log
   f <- input$f_log
   D <- input$d_log
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- (1-omega)/((1-omega) + omega*time^gamma)

   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="probability", main="Survival Function")
   plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="Measurement Occasion", ylab="Probability", main="Hazard Function")
 }, height=600)
 
 ### Gompertz function
 output$gompplots <- renderPlot({
   
   validate(need(input$d_gomp>=3,"Duration of the study should be at least three."))
   if(input$omega_gomp < 0 || input$omega_gomp > 1) {
     validate("`Omega` should be between zero and one.")
   }
   validate(need(input$gamma_gomp != 0, "`Gamma` cannot be zero"))
   
   omega <- input$omega_gomp
   gamma <- input$gamma_gomp
   f <- input$f_gomp
   D <- input$d_gomp
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- exp((log(1-omega)/(exp(gamma)-1))*(exp(gamma*time)-1))

   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="probability", main="Survival Function")
   plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="Measurement Occasion", ylab="Probability", main="Hazard Function")
 }, height=600)
 
### Non-parametric
 
 output$hand.plots <- renderPlot({
   if(input$option == "opt1") {
       t.points.hand <- as.numeric(unlist(strsplit(input$meas.occ,",")))
       surv.hand <- as.numeric(unlist(strsplit(input$survival.hand,",")))
       validate(need(length(t.points.hand) == length(surv.hand),"The number of measurement occasions must be equal to the length of S(j)."))
       validate(need(is.unsorted(rev(surv.hand)) == FALSE, "The expected proportion of participants remaining can never increase as we assume that participants who dropped out cannot re-enter the study at a later point in time."))
       n <- length(t.points.hand)
       hazard <- rep(NA, n)
       for(i in 1:(n-1)) {
         hazard[i+1] <- round((surv.hand[i]-surv.hand[i+1])/surv.hand[i], digits=7)
       }
       par(mfrow=c(2,1))
       plot(t.points.hand, surv.hand, type="b", xaxt="n", pch=16, ylim=c(0,1), xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Probability", main="Your own Survival Function")
       axis(side=1, at=t.points.hand)
       plot(t.points.hand, hazard, type="b", xaxt="n", pch=16, xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Probability", main="The Resulting Hazard Function")
       axis(side=1, at=t.points.hand)
       
   } else if(input$option == "opt2") {
       t.points.hand <- as.numeric(unlist(strsplit(input$meas.occ,",")))
       haz.hand <- as.numeric(unlist(strsplit(input$hazard.hand,",")))
       validate(need(length(t.points.hand) == length(haz.hand)+1,"The number of measurement occasions must be equal to the length of h(j) + 1."))
       n <- length(t.points.hand)
       survival <- c(1, rep(NA, (n-1)))
       for(i in 2:n){
         survival[i] <- (1-haz.hand[i-1]) * survival[i-1]
       }
       par(mfrow=c(2,1))
       plot(t.points.hand, survival, type="b", xaxt="n", pch=16, ylim=c(0,1), xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Probability", main="The Resulting Survival Function")
       axis(side=1, at=t.points.hand)
       plot(t.points.hand, c(NA, haz.hand), type="b", xaxt="n", pch=16, xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="probability", main="Your own Hazard Function")
       axis(side=1, at=t.points.hand)
   }
 }, height=600)
   
}
