# Server part of the shiny app

server <- function(input, output) {

### Welcome page
output$welcome_plot1 <- renderPlot({
  
  x <- seq(0, 10)
  time <- seq(0, 1, length=length(x))
  omega <- .8
  gamma <- .7
  y <- (1-omega)^time^gamma
  
  plot(x, y,
       type = "b",
       lwd = 2, pch=21, bg="#ffd800", cex=1.5, 
       ylim = c(0, 1),
       xlab = "Measurement occasion",
       ylab = "Survival probability",
       main = "Survival modelled with the Weibull function")
  
})

output$welcome_plot2 <- renderPlot({
  
  x <- seq(0, 10)
  time <- seq(0, 1, length=length(x))
  omega <- .8
  gamma <- .7
  remain <- (1-omega)^time^gamma
  
  hazard <- rep(NA, length(x))
  for(i in 1:(length(x)-1)) {
    hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
  }
  
  plot(x, hazard,
       type = "b",
       lwd = 2, pch=21, bg="#ffd800", cex=1.5, 
       ylim = c(0, 1),
       xlab = "Measurement occasion",
       ylab = "Hazard probability",
       main = "Example: Decreasing hazard over time")
  
})
  
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
      hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
    }
    
  par(mfrow=c(2,1))
  plot(t.points,remain,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="Survival probability", main="Survival Function")
  grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
  plot(t.points,hazard,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,D), xlab="Measurement Occasion", ylab="Hazard probability", main="Hazard Function")
  grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
  })
 
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
     hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="Survival probability", main="Survival Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
   plot(t.points,hazard,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,D), xlab="Measurement Occasion", ylab="Hazard probability", main="Hazard Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
})
 
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
     hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="Survival probability", main="Survival Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
   plot(t.points,hazard,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,D), xlab="Measurement Occasion", ylab="Hazard probability", main="Hazard Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
 })
 
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
     hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="Survival probability", main="Survival Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
   plot(t.points,hazard,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,D), xlab="Measurement Occasion", ylab="Hazard probability", main="Hazard Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
 })
 
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
     hazard[i] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,D), xlab="Measurement Occasion", ylab="Survival probability", main="Survival Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
   plot(t.points,hazard,type="b", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,D), xlab="Measurement Occasion", ylab="Hazard probability", main="Hazard Function")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
 })
 
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
         hazard[i] <- round((surv.hand[i]-surv.hand[i+1])/surv.hand[i], digits=7)
       }
       par(mfrow=c(2,1))
       plot(t.points.hand, surv.hand, type="b", xaxt="n", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Survival probability", main="Your own Survival Function")
       grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
       axis(side=1, at=t.points.hand)
       plot(t.points.hand, hazard, type="b", xaxt="n", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Hazard probability", main="The Resulting Hazard Function")
       grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
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
       plot(t.points.hand, survival, type="b", xaxt="n", lwd=2, pch=21, bg="#ffd800", cex=1.5, ylim=c(0,1), xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Survival probability", main="The Resulting Survival Function")
       grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
       axis(side=1, at=t.points.hand)
       plot(t.points.hand, c(NA, haz.hand), type="b", xaxt="n", lwd=2, pch=21, bg="#ffd800", cex=1.5, xlim=c(0,max(t.points.hand)), xlab="Measurement Occasion", ylab="Hazard probability", main="Your own Hazard Function")
       grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
       axis(side=1, at=t.points.hand)
   }
 })
 
### Click Plot
   
 values <- reactiveValues(
   x = 0,
   y = 1
 )
 
 # input values by clicking
 observeEvent(input$plot_click, {
   if (input$plot_click$x > input$x_max) {
     showNotification("Click is outside x-axis range", type = "error")
     return()
   }
   if (length(values$x) >= 50) {
     showNotification("Maximum of 50 points reached", type = "warning")
     return()
   }
   
   new_x <- round(input$plot_click$x)
   new_y <- abs(input$plot_click$y)
   
   # prevent negative x
   if (new_x < 0) {
     showNotification("x must be ≥ 0", type = "error")
     return()
   }
   
   tol <- 1e-6
   # SPECIAL CASE: x ≈ 0
   if (abs(new_x) < tol) {
     idx <- which(abs(values$x) < tol)
     if (length(idx) > 0) {
       values$y[idx] <- 1
     } else {
       values$x <- c(values$x, 0)
       values$y <- c(values$y, 1)
     }
     return()
   }
   
   # normal duplicate handling
   idx <- which(abs(values$x - new_x) < tol)
   if (length(idx) > 0) {
     values$y[idx] <- min(values$y[idx], new_y)
   } else {
     values$x <- c(values$x, new_x)
     values$y <- c(values$y, new_y)
   }
 })
 
 # reset button
 observeEvent(input$reset_plot, {
   values$x <- 0
   values$y <- 1
 })
 
 # if x_max is changed
 observeEvent(input$x_max, {
   keep <- values$x <= input$x_max
   values$x <- values$x[keep]
   values$y <- values$y[keep]
 })
 
 # plot output
 output$click_plot <- renderPlot({
   
   x_max <- input$x_max
   
   plot(NA, NA,
        xlim = c(0, x_max),
        ylim = c(0, 1),
        type = "b",
        xlab = "Measurement occasion",
        ylab = "Survival probability")
   grid(nx = NULL, ny = NULL, col = "gray90", lty = 1)
   
   if (length(values$x) > 0) {
     
     ord <- order(values$x)
     x <- values$x[ord]
     y <- cummin(values$y[ord])
     
     points(x, y,
            pch = 21,
            bg = "#ffd800",   # fill color
            col = "black",    # border color
            cex = 2,
            lwd = 1.5)
     
     lines(x, y, lwd = 2)
   }
 })
 
 # Table of clicked values
 output$click_table <- renderTable({
   
   if (length(values$x) == 0) return(NULL)
   
   ord <- order(values$x)
   x <- values$x[ord]
   y <- cummin(values$y[ord])
   
   x_display <- formatC(x, format = "fg", digits = 6)
   y_display <- round(y*100, 1)
   
   df <- rbind(
     "timepoints" = x_display,
     "% of sample remaining" = y_display
   )
   
   colnames(df) <- paste0("occasion ", seq_along(x))
   
   df
 }, rownames = TRUE)
 
}
