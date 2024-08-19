library(shiny)
library(shinydashboard)
library(shinythemes)
library(DT)
library(bslib) 

# Here, we will assume that individuals who dropped out do not re-enter the study at a later point in time. Note we assume that the chance of drop-out can depend on the amount of time elapsed but not on the number of measurement occasions.

ui <- page_navbar(title ="Attrition", bg = "#ffd800",
                  nav_panel(title="Introduction",
                            withMathJax(),

                            tags$head(tags$script(type = "text/javascript", src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js")),
                            tags$h1("Modelling attrition patterns in longitudinal trials using survival functions"),
                            tags$p("Attrition refers to the drop-out of participants during the course of a study. This results in a loss of information that can have large negative effects on statistical power. Therefore, one should take expected drop-out rates into account when performing sample size determination (SSD) to avoid underpowered studies (Moerbeek, 2020). For SSD for longitudinal trials, this can be done by modelling the expected attrition pattern via a survival function. To do this, we make a number of assumptions:"),
                            tags$ol(
                              tags$li("Individuals who dropped out do not re-enter the study at a later point in time."),
                              tags$li("Although the underlying time variable is continuous, we can only observe dropout when a measurement is taken. This means that the time variable in the plots is discrete."),
                              tags$li("The chance of drop-out depends on the amount of time elapsed since the start of the study but not on the number of measurement occasions taken so far.")
                            ),
                            tags$h4("Survival probability functions and hazard probability functions"),
                            tags$p("The survival function denotes the probability of remaining in the study up until measurement occasion j and is denoted as $$S(j)=P(J \\geq j),$$ where J represents the number of measurement occasions that an individual has completed without dropping out. The graph of this function can be interpreted as the expected proportion of individuals who dropped out at measurement occasion j. The hazard function h(j), on the other hand, is a conditional discrete distribution representing the probability of an individual dropping out at measurement occasion j, given that they have remained in the study up until j-1. The graph of this function can be interpreted as the probability of an individual dropping out between measurement occasion j-1 and j."),
                            tags$h4("Exponential survival function"),
                            tags$p("After rescaling time to represent the proportion of time elapsed [0,1] in a study and replacing the parameter lambda with \\–log(1-omega)\\ for better interpretability, the exponential survival function is defined as 
                                   \\S(t)=(1-omega)^t\\. Omega [0,1] represents the proportion of participants who drop out at some point during the course of the study. When modelling attrition with the exponential function, one assumes that the hazard rate is constant, meaning that the chances of dropping out are the same at each point in time. This is illustrated by the fact that the slope of the hazard function always equals zero, \\h’(t)=0\\."),
                            tags$h4("Weibull survival function"),
                            tags$p("The Weibull survival function can be seen as a generalization of the exponential function. It is defined as \\(1-omega)^t^gamma\\, where gamma [0,INF] denotes the slope of the hazard function. If \\gamma=0\\, then the Weibull function reduces to the exponential function. For \\gamma>1\\, the slope of the hazard function is positive, meaning that attrition is higher towards the end of the study. For \\gamma<1\\, the slope of the hazard function is negative, meaning that drop-out is higher at the beginning of the study."),
                            tags$h4("This ShinyApp"),
                            tags$p("At the top of the page, you can navigate to the survival function suitable to model the expected pattern of attrition in you trial. The plots visualize both the survival and the hazard function, giving a visual impression of the effect of altering the parameters of the functions. The function to perform Bayesian SSD for longitudinal trials can be found in this GitHub repository: https://github.com/ulrichlosener/BayesianSSD_attrition."),
                            tags$h4("References"),
                            tags$p("Moerbeek, M. (2020). The cluster randomized crossover trial: The effects of attrition in the AB/BA design and how to account for it in sample size calculations. Clinical Trials, 17(4), 420-429.")
                            ),
                  nav_panel(title="Weibull", 
                            layout_sidebar(sidebar = sidebar(
                              numericInput("D1", "Duration of study", 5),
                              numericInput("f1", "Frequency of observation", 1),
                              numericInput("omega1", "Omega", .5, min = 0, max = 1, step = .1),
                              numericInput("gamma", "Gamma", 3, min = 0, max = 100, step = .5)
                            ),
                            fluidRow(column(width=6, plotOutput(outputId = "weibullplots")),
                                     column(width=6,
                                            tags$div(
                                              tags$h1("Explanation of parameters"),
                                              tags$h4("Duration of study"),
                                                tags$ul(
                                                  tags$li("The total duration of the study in days, weeks, months, etc."),
                                                  tags$li("Must be at least three in this case")
                                                ),
                                              tags$h4("Frequency of observation"),
                                                tags$ul(
                                                  tags$li("The number of observations taken per unit in time"),
                                                  tags$li("If measurements are taken every two weeks, frequency = 0.5; if measurements are taken twice a day, frequency = 2"),
                                                  tags$li("Note that here, we assume the measurements to be equidistant, i.e., equally spaced in time ")
                                                ),
                                              tags$h4("Omega"),
                                                tags$ul(
                                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                                ),
                                              tags$h4("Gamma"),
                                                tags$ul(
                                                  tags$li("Gamma: Hazard rate, determines whether dropout in concentrated towards the beginning (gamma < 1) or the end (gamma > 1) of a study"),
                                                  tags$li("If gamma = 1, then the Weibull function equals the exponential function")
                                                )
                                            )
                                      )
                            )
                            )
                  ),
                  
                  nav_panel(title="Exponential",
                            layout_sidebar(sidebar = sidebar(
                              numericInput("D2", "Duration of study", 5),
                              numericInput("f2", "Frequency of observation", 1),
                              numericInput("omega2", "Omega", .5, min = 0, max = 1, step = .1)
                            ),
                            fluidRow(column(width=6, plotOutput(outputId = "exponentialplots")),
                                     column(width=6, 
                                            tags$div(
                                              tags$h1("Explanation of parameters"),
                                              tags$h4("Duration of study"),
                                                tags$ul(
                                                  tags$li("The total duration of the study in days, weeks, months, etc."),
                                                  tags$li("Must be at least three in this case")
                                                ),
                                              tags$h4("Frequency of observation"),
                                                tags$ul(
                                                  tags$li("The number of observations taken per unit in time"),
                                                  tags$li("If measurements are taken every two weeks, frequency = 0.5; if measurements are taken twice a day, frequency = 2"),
                                                  tags$li("Note that here, we assume the measurements to be equidistant, i.e., equally spaced in time ")
                                                ),
                                              tags$h4("Omega"),
                                                tags$ul(
                                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                                )
                                            )
                                    )
                            )
                            )
                  ),
                  
                  nav_panel(title="More survival functions", p("Coming soon!")),
                  nav_spacer()
)








