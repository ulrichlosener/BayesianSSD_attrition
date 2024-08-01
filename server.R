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
  plot(t.points,remain,type="b",pch=16,ylim=c(0,1),xlim=c(0,D),xlab="time",ylab="probability",main="Survival function")
  plot(t.points,hazard,type="b",pch=16,xlim=c(0,D),xlab="time",ylab="probability", main="Hazard function")
    }, height=600)
 
 output$weibullinfo <- renderText({
   "The Weibull distribution is a generalization of the Exponential distribution in a sense that the hazard rate is not constrained to be constant over time.
   The parameter `omega` represents the proportion of individuals who drop out of the trial at some point.
   The parameter `gamma` represents the linear change in dropout propbability over time. 
   If gamma < 1, individuals are more likely to drop out at the beginning of the trial.
   If gamma > 1, individuals are more likely to drop out at the end of the trial.
   If gamma = 1, the dropout risk is constant (in this case, the Weibull distribution is equivalent to the exponential distribution)
   "
 })
 
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
   plot(t.points,remain,type="b",pch=16,ylim=c(0,1),xlim=c(0,D),xlab="time",ylab="probability",main="Survival function")
   plot(t.points,hazard,type="b",pch=16,xlim=c(0,D),xlab="time",ylab="probability", main="Hazard function")
 }, height=600)
 
 output$exponentialinfo <- renderText({
   "The exponential function is a special case of the Weibull function where the hazard rate is constant, i.e., gamma=1. The parameter `omega` represents the proportion of individuals who drop out of the trial at some point."
 })
}
