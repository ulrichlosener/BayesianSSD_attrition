server <- function(input, output) {

 output$weibullplots <- renderPlot({
    
    validate(need(input$D1>=3,"Duration of the study should be at least three."))
    if(input$omega1<0 || input$omega1>1) {
      validate("`Omega` should be between zero and one.")
    }
    validate(need(input$gamma>=0,"`Gamma` should be positive."))
    
    omega <- input$omega1
    gamma <- input$gamma
    f <- input$f1
    D <- input$D1
    t.points <- seq(0, D, by=1/f)
    time <- seq(0, 1, length=length(t.points))
    remain <- (1-omega)^time^gamma
    n <- length(remain)
    
    hazard <- rep(NA, n)
    for(i in 1:(n-1)) {
      hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
    }
    
  par(mfrow=c(2,1))
  plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="measurement occasion", ylab="probability", main="Survival function")
  plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="measurement occasion", ylab="probability", main="Hazard function")
    }, height=600)
 
 output$exponentialplots <- renderPlot({
   validate(need(input$D2>=3,"Duration of the study should be at least three."))
   if(input$omega2<0 || input$omega2>1) {
     validate("`Omega` should be between zero and one.")
   }
   
   omega <- input$omega2
   f <- input$f2
   D <- input$D2
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- (1-omega)^time
   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b", pch=16, ylim=c(0,1), xlim=c(0,D), xlab="measurement occasion", ylab="probability", main="Survival function")
   plot(t.points,hazard,type="b", pch=16, xlim=c(0,D), xlab="measurement occasion", ylab="probability", main="Hazard function")
 }, height=600)
 
 output$hand.plots <- renderPlot({
   if(input$option == "opt1") {
       t.points.hand <- as.numeric(unlist(strsplit(input$meas.occ,",")))
       surv.hand <- as.numeric(unlist(strsplit(input$survival.hand,",")))
       n <- length(t.points.hand)
       hazard <- rep(NA, n)
       for(i in 1:(n-1)) {
         hazard[i+1] <- round((surv.hand[i]-surv.hand[i+1])/surv.hand[i], digits=7)
       }
       par(mfrow=c(2,1))
       plot(t.points.hand, surv.hand, type="b", xaxt="n", pch=16, ylim=c(0,1), xlim=c(0,max(t.points.hand)), xlab="measurement occasion", ylab="probability", main="Your own Survival function")
       axis(side=1, at=t.points.hand)
       plot(t.points.hand, hazard, type="b", xaxt="n", pch=16, xlim=c(0,max(t.points.hand)), xlab="measurement occasion", ylab="probability", main="Your own Hazard function")
       axis(side=1, at=t.points.hand)
   } else if(input$option == "opt2") {
     plot(c(1,2,3), c(1,2,3))
   }
 }, height=600)
   
}
