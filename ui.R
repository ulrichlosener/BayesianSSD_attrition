# user interface of the shiny app

library(shiny)
library(shinydashboard)
library(shinythemes)
library(DT)
library(bslib) 


ui <- page_navbar(title ="Attrition", bg = "#ffd800",
                  nav_panel(title="Introduction",
                            tags$head(tags$script(type = "text/javascript", src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js")),
                            tags$h1("Modelling attrition patterns in longitudinal trials using survival functions"),
                            tags$p("Attrition refers to the drop-out of participants during the course of a study. This results in a loss of information that can have large negative effects on statistical power. Therefore, one should take expected drop-out rates into account when performing sample size determination (SSD) to avoid underpowered studies (Moerbeek, 2020). For SSD for longitudinal trials, this can be done by modelling the expected attrition pattern via a survival function. To do this, we make a number of assumptions:"),
                            tags$ol(
                              tags$li("Individuals who dropped out do not re-enter the study at a later point in time."),
                              tags$li("Although the underlying time variable is continuous, we can only observe dropout when a measurement is taken. This means that the time variable in the plots is discrete."),
                              tags$li("The chance of drop-out depends on the amount of time elapsed since the start of the study but not on the number of measurement occasions taken so far.")
                            ),
                            tags$h4("Survival probability functions and hazard probability functions"),
                            tags$p("The survival function denotes the probability of remaining in the study up until measurement occasion j and is denoted as $$S(j)=P(J \\geq j),$$ where J represents the number of measurement occasions that an individual has completed without dropping out. Thus, the survival probability function is a discrete distribution showing the probability of having remained in the study up until measurement occasion j (where j = 0, 1, 2, ..., n). The graph of this function can be interpreted as the expected proportion of individuals who dropped so far for each measurement occasion j. The hazard function, on the other hand, is a conditional discrete distribution representing the probability of an individual dropping out at measurement occasion j, given that they have remained in the study up until j - 1 and is defined as $$h(j)=\\frac{S(j)-S(j+1)}{S(j)}.$$ The graph of this function can therefore be interpreted as the probability of an individual dropping out between measurement occasion j - 1 and j."),
                            tags$h4("Exponential survival function"),
                            tags$p("After rescaling j to represent the proportion of measurement occasions taken in a study, the exponential survival function is defined as $$S_{exponential}(j)=(1-\\omega)^j,$$ where omega [0,1] represents the proportion of participants who remain in the study without dropping out. When modelling attrition with the exponential probability function, one assumes that the hazard probability is constant, meaning that the chances of having dropped out are the same at each measurement occasion j. This is illustrated by the fact that the slope of the exponential hazard probability function is always equals zero, h’(j) = 0."),
                            tags$h4("Weibull survival function"),
                            tags$p("The Weibull survival function can be seen as a generalization of the exponential function where j is given an exponent, gamma. It is defined as $$S_{weibull}(j)={(1-\\omega)^j}^\\gamma$$ where gamma [0,INF] denotes the slope of the hazard function. It can be seen that if gamma = 0, then the Weibull function reduces to the exponential function. For gamma > 1, the slope of the hazard function is positive, meaning that attrition is higher towards the end of the study. For gamma < 1, the slope of the hazard function is negative, meaning that drop-out is higher at the beginning of the study."),
                            tags$h4("This ShinyApp"),
                            tags$p("At the top of the page, you can navigate to the survival function suitable to model the expected pattern of attrition in you trial. The plots visualize both the survival and the hazard function, giving a visual impression of the effect of altering the parameters of the functions. The function to perform Bayesian SSD for longitudinal trials can be found in this GitHub repository: https://github.com/ulrichlosener/BayesianSSD_attrition."),
                            tags$h4("References"),
                            tags$p("Moerbeek, M. (2020). The cluster randomized crossover trial: The effects of attrition in the AB/BA design and how to account for it in sample size calculations. Clinical Trials, 17(4), 420-429."),
                            tags$br("")
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
                  nav_panel(title="Non-parametric Survival function", 
                            layout_sidebar(sidebar = sidebar(
                              radioButtons("option", "Choose an option:", 
                                           choices = list("Enter the survival probability function, i.e., the expected proportion of participants remaining in the study for each measurement occasion" = "opt1", 
                                                          "Enter the hazard probability function, i.e., the probability of dropping out between j and j-1" = "opt2"), selected = "opt1"),
                              conditionalPanel(condition = "input.option == 'opt1'", 
                                               textInput("meas.occ", "Measurement occasions", "0, 1, 2, 5, 8, 9, 10"),
                                               textInput("survival.hand", "Expected proportion of participants remaining in the study at each measurement occasion, S(j)", "1, .95, .9, .8, .7, .65, .6")
                                               ),
                              conditionalPanel(condition = "input.option == 'opt2'",
                                               textInput("meas.occ", "Measurement occasions (comma delimited)", "0, 1, 2, 5, 8, 9, 10"),
                                               textInput("hazard.hand", "Expected risk of dropping between measurement occasion j and j-1, h(j)", ".1, .15, .2, .25, .3, .5"),
                                               )
                            ), 
                            fluidRow(column(width=6, plotOutput(outputId = "hand.plots")),
                                     column(width=6, 
                                            tags$h4("Choosing an option"),
                                            tags$p("Here, you can choose between specifying your own survival probability function (option 1) or your own hazard probability function (option 2). Whatever is more intuitive for you is the better option. Either way, you will see both the plot for the survival probability function as well as the plot for the hazard probability function."),
                                            tags$h4("Option 1: Specifying a survival probability function"),
                                            tags$p("First, the measurement occasions need to be specified. These can be non-equidistant, meaning that the time interval between observations does not need to be constant. If you measure, say in week 0 (baseline), in week 1, in week 5, and in week 10, then you should input '0, 1, 5, 10' for the measurement occasions. Next, you input the expected proportion of participants still remaining in the study for each of these measurement occasions (starting with 1 for the first measurement). For example, if you expect 10% of participants to drop out after the first measurement, 15% after the second measurement, and another 30% after the third measurement, you should input '1, .9, .85, .7'. Note that the first value will always be 1 as there can be no dropout before the first measurement occasion."),
                                            tags$h4("Option 2: Specifying a hazard probability function"),
                                            tags$p("The specification of the measurement occasions is the same as for option 1. Next, you input the expected risk of dropout between measurement occasion j and j - 1 for each j = 1, 2, ..., n. Note that the value for occasion j = 0 is undefined because there can be no dropout before the first measurement occasion. Using the same example measurement occasions as in option 1, if you expect the risk of dropping out in the first week to be .10, between week 1 and week 5 to be .20, and between week 5 and week 10 to be .40, your input should be '.1, .2, .4'. Note that the number of elements in expected risk need to be one less that the amount of measurement occasions as h(j=0) is undefined.")
                                     ))
                                  )               
                            ),
                  nav_panel(title="More?", 
                            tags$h4("If you are missing information or a certain functionality in this ShinyApp, please feel free to contact me!"),
                            tags$p("Via e-mail:"),
                            tags$p("u.c.losener1@uu.nl"),
                            tags$p("or see the following link for more information."),
                            tags$a(href="https://www.uu.nl/staff/UCLosener2", "Ulrich Lösener"),
                            ),
                  nav_spacer(),
                  nav_menu(
                    title = "Links",
                    nav_item(tags$a("GitHub", href = "https://github.com/ulrichlosener/BayesianSSD_attrition")),
                    nav_item(tags$a("Contact", href = "https://www.uu.nl/staff/UCLosener2/Contact"))
                  )
  )








