server <- function(input, output) {

 output$weibullplots <- renderPlot({
    
    validate(need(input$D>=3,"Duration of the study should be at least three."))
    if(input$omega<0 || input$omega>1) {
      validate("`Omega` should be between zero and one.")
    }
    validate(need(input$gamma>=0,"`Gamma` should be positive."))
    
    omega <- input$omega
    gamma <- input$gamma
    f <- input$f
    D <- input$D
    t.points <- seq(0, D, by=1/f)
    time <- seq(0, 1, length=length(t.points))
    remain <- (1-omega)^time^gamma
    n <- length(remain)
    
    hazard <- rep(NA, n)
    for(i in 1:(n-1)) {
      hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
    }
    
  par(mfrow=c(2,1))
  plot(t.points,remain,type="b",pch=16,ylim=c(0,1),xlim=c(0,D),xlab="time",ylab="probability",main="Survival function")
  plot(t.points,hazard,type="b",pch=16,xlim=c(0,D),xlab="time",ylab="probability", main="Hazard function")
    })
 
 output$weibullinfo <- renderText({
   "Here some very informative stuff about the Weibull function"
 })
 
 output$exponentialplots <- renderPlot({
   validate(need(input$D>=3,"Duration of the study should be at least three."))
   if(input$omega<0 || input$omega>1) {
     validate("`Omega` should be between zero and one.")
   }
   
   omega <- input$omega
   f <- input$f
   D <- input$D
   t.points <- seq(0, D, by=1/f)
   time <- seq(0, 1, length=length(t.points))
   remain <- (1-omega)^time
   n <- length(remain)
   
   hazard <- rep(NA, n)
   for(i in 1:(n-1)) {
     hazard[i+1] <- round((remain[i]-remain[i+1])/remain[i], digits=7)
   }
   
   par(mfrow=c(2,1))
   plot(t.points,remain,type="b",pch=16,ylim=c(0,1),xlim=c(0,D),xlab="time",ylab="probability",main="Survival function")
   plot(t.points,hazard,type="b",pch=16,xlim=c(0,D),xlab="time",ylab="probability", main="Hazard function")
    })
 
   output$exponentialinfo <- renderText({
     "Here some very informative stuff about the exponential function"
   })
}
