# user interface of the shiny app

library(shiny)
library(shinydashboard)
library(shinythemes)
library(DT)
library(bslib) 

ui <- page_navbar(title ="Attrition", bg = "#ffd800", theme = bslib::bs_theme(version = 5, base_font = bslib::font_google("Roboto")),
                  nav_panel(title = "Welcome",
                    div(
                      class = "container-fluid",
                      style = "max-width: 900px; margin: auto;",
                      # Title
                      tags$h1("Welcome to the Attrition App!"),
                      tags$p(
                        "This application helps you explore and model attrition patterns ",
                        "in longitudinal studies using survival functions."
                      ),
                      
                      tags$br(),
                      
                      # Section 1
                      div(
                        class = "card p-3 shadow-sm",
                        tags$h4("What can you do here?"),
                        tags$ul(
                          tags$li("Explore different parametric and non-parametric survival models"),
                          tags$li("Adjust parameters and immediately see their effect on the curve"),
                          tags$li("Draw your own survival curve manually"),
                          tags$li("See how survival and hazard curves are connected")
                        )
                      ),
                      
                      tags$br(),
                      
                      # Section 2
                      div(
                        class = "card p-3 shadow-sm",
                        tags$h4("How to get started"),
                        tags$ol(
                          tags$li("Choose a model from the navigation bar"),
                          tags$li("Adjust parameters using the input fields"),
                          tags$li("Inspect how the survival and hazard functions change"),
                          tags$li("Or go to 'Survival by hand' to draw your own curve")
                        )
                      ),
                      
                      tags$br(),
                      
                      # Example plot
                      tags$h4("Example: Weibull survival function"),
                      tags$ul(
                        tags$li("Shows the probability of 'surviving' in the sample up to a certain point"),
                        tags$li("80% of participants drop out at some point (omega = .8)"),
                        tags$li("Chance of dropout is higher at the beginning of the study (gamma = 0.7)")
                      ),
                      div(
                        style = "height:40vh;",
                        plotOutput("welcome_plot1", height = "100%")
                      ),
                      tags$h4("The corresponding hazard function"),
                      tags$ul(
                        tags$li("Shows the dropout risk between the current and the subsequent measurement"),
                        tags$li("In this example, dropout risk decreases over time (because we set gamma < 1)")
                      ),
                      div(
                        style = "height:40vh;",
                        plotOutput("welcome_plot2", height = "100%")
                      )
                    )
                  ),
                  nav_panel(title="Background",
                                tags$head(tags$style(HTML("
                                  @media (max-width: 768px) {
                                    h1 { font-size: 20px; }
                                    h4 { font-size: 16px; }
                                    p, li { font-size: 14px; }
                                  }
                                  .shiny-input-container {
                                    width: 100% !important;
                                  }
                                ")),
                                tags$script(type = "text/javascript", src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js")),
                                tags$h1("Modelling attrition patterns in longitudinal trials using survival functions"),
                                p("Attrition refers to the drop-out of participants during the course of a study. This results in a loss of information that can have large negative effects on statistical power. Therefore, one should take expected drop-out rates into account when performing sample size determination (SSD) to avoid underpowered studies (Moerbeek, 2020). Suppose we want to measure N individuals n times for the purpose of testing hypotheses about a potential treatment effect over time and want to take expected attrition into account in our SSD. This can be done by modelling the expected attrition pattern via a survival function. To do this, we make a number of assumptions:"),
                                tags$ol(
                                  tags$li("Individuals who dropped out do not re-enter the study at a later point in time."),
                                  tags$li("Although the underlying time variable is continuous, we can only observe dropout when a measurement is taken. This means that the time variable in the plots is discrete."),
                                  tags$li("The chance of drop-out depends on the amount of time elapsed since the start of the study but not on the number of measurement occasions taken so far.")
                                ),
                                tags$h4("Survival probability functions and hazard probability functions"),
                                p("The discrete survival function denotes the probability of remaining in the study up until measurement occasion j and is denoted as $$S(j)=P(J \\geq j),$$ where J represents the number of measurement occasions that an individual has completed without dropping out. Thus, the survival probability function is a discrete distribution showing the probability of having remained in the study up until measurement occasion j (where j = 0, 1, 2, ..., n). The graph of this function can be interpreted as the expected proportion of individuals who dropped so far for each measurement occasion j. The hazard function, on the other hand, is a conditional discrete distribution representing the probability of an individual dropping out at measurement occasion j, given that they have remained in the study up until j - 1 and is defined as $$h(j)=\\frac{S(j)-S(j+1)}{S(j)}.$$ The graph of this function can therefore be interpreted as the probability of an individual which has thus far remained in the study dropping out between measurement occasion j - 1 and j. Note that for convenience, we rescale j to represent the proportion of measurement occasions taken in a study from this point forward."),
                                tags$h4("Weibull survival function"),
                                p("The Weibull survival function is defined as $$S_{weibull}(j)={(1-\\omega)^j}^\\gamma$$ where omega [0, 1] represents the proportion of participants who drop out of the study and gamma [0, inf] denotes the slope of the hazard function. If gamma = 0, the Weibull function reduces to the exponential function and has a constant hazard rate. For gamma > 1, the slope of the hazard function is positive, meaning that attrition is higher towards the end of the study. For gamma < 1, the slope of the hazard function is negative, meaning that drop-out is higher at the beginning of the study."),
                                tags$h4("Modified Weibull survival function"),
                                p("The Modified Weibull survival function is defined as $$S_{weibull}(j)=e^{j \\gamma e^{\\kappa (j-1)}log(1 - \\omega)}$$ where omega [0, 1] represents the proportion of participants who drop out of the study and gamma and kappa [0, inf] are shape parameters. If gamma > 0, the hazard monotonically increases, meaning that attrition is higher towards the end of the study. If gamma < 0, the hazard rate has a bathtub shape with turning point at $$j_{min}=(\\sqrt\\gamma - \\gamma)/ \\kappa.$$ meaning that dropout is highest in the beginning and at the end of a study."),
                                tags$h4("Log-Logistic survival function"),
                                p("The Log-Logistic survival function is defined as $$S_{log-logistic}(j)=\\frac{1 - \\omega}{1 - \\omega + \\omega j^\\gamma}$$ where omega [0, 1] represents the proportion of participants who drop out of the study and gamma [0, inf] denotes the shape of the hazard function. If gamma < 1, the hazard rate monotonically decreases, meaning that attrition is concentrated at towards the beginning of the study. For gamma < 1, the hazard rate is bell-shaped with its peak at $$j_{max}=\\frac{(\\gamma - 1)^{1 / \\gamma}}{(\\omega / (1 - \\omega))^{1/\\gamma}}$$ meaning that dropout is concentrated somewhere halfway through the study."),
                                tags$h4("Linear-Exponential survival function"),
                                p("The Linear-Exponential survival function is defined as $$S_{linear-exponential}(j)=e^{(0.5 \\gamma + log(1 - \\omega))j - 0.5 \\gamma j^2}$$ where omega [0, 1] represents the proportion of participants who drop out of the study and gamma [-1, 1] determines the shape of the distribution. If gamma > 0, the slope of the hazard function is linearly positive, meaning that attrition is higher towards the end of the study. If gamma < 0, the slope of the hazard function is linearly negative, meaning that drop-out is higher at the beginning of the study."),
                                tags$h4("Gompertz survival function"),
                                p("The Gompertz survival function is defined as $$S_{gompertz}(j)=e^{(log(1 - \\omega)/(e^\\gamma - 1))(e^{\\gamma j} - 1)}$$ where omega [0, 1] represents the proportion of participants who drop out of the study and gamma ≠ 0 denotes the shape of the hazard function. If gamma > 0, the slope of the hazard function is exponentially positive, meaning that attrition is higher towards the end of the study. If gamma < 0, the slope of the hazard function is exponentially negative, meaning that drop-out is higher at the beginning of the study."),
                                tags$h4("This ShinyApp"),
                                p("At the top of the page, you can navigate to the survival function suitable to model the expected pattern of attrition in you trial. The plots visualize both the survival and the hazard function, giving a visual impression of the effect of altering the parameters of the functions. The open-access R function to perform Bayesian SSD for longitudinal trials with attrition can be found in this GitHub repository: https://github.com/ulrichlosener/BayesianSSD_attrition."),
                                tags$h4("References"),
                                p("Moerbeek, M. (2020). The cluster randomized crossover trial: The effects of attrition in the AB/BA design and how to account for it in sample size calculations. Clinical Trials, 17(4), 420-429."),
                                tags$br("")
                            ),
                  nav_panel(title="Survival by hand",
                            tags$h4("Here, you can manually draw the points of your own non-parametric survival function"),
                            tags$p("Click inside the plot to add points, the exact values are displayed in the table at the bottom."),
                            numericInput("x_max", "Adjust the study duration (x-axis maximum)", 
                                         value = 100, min = 1, step = 1),
                            div(
                              style = "height:75vh;",  
                              plotOutput("click_plot", click = "plot_click", height = "100%")
                            ),
                            actionButton("reset_plot", "Reset points", style="color: black; background-color: #ffd800"),
                            div(
                              style = "overflow-x: auto;",
                              tableOutput("click_table")
                            )
                  ),
                  nav_panel(title="Weibull", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always",  
                              numericInput("d_weib", "Duration of study", 10),
                              numericInput("f_weib", "Frequency of observation", 2),
                              numericInput("omega_weib", "Omega", .5, min = 0, max = 1, step = .1),
                              numericInput("gamma_weib", "Gamma", 3, min = 0, max = 100, step = .5)
                            ),
                            fluidRow(
                              column(width = 12, class = "col-md-6", 
                                     div(
                                       class = "d-flex flex-column",
                                       style = "height:90vh;",
                                       plotOutput("weibullplots", height = "100%")
                                     )
                              ),
                              column(width = 12, class = "col-md-6", 
                                     tags$div(tags$div(style = "padding:15px;",
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
                                                  tags$li("[0, 1]"),
                                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                                ),
                                                tags$h4("Gamma"),
                                                tags$ul(
                                                  tags$li("[0, inf]"),
                                                  tags$li("Hazard rate, determines whether dropout in concentrated towards the beginning (gamma < 1) or the end (gamma > 1) of a study"),
                                                  tags$li("If gamma = 1, then the Weibull function equals the exponential function")
                                                )
                                              )
                                     )
                                     )
                            )
                            )
                  ),
                  nav_panel(title="Modified-Weibull", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always", 
                              numericInput("d_mod_weib", "Duration of study", 10),
                              numericInput("f_mod_weib", "Frequency of observation", 2),
                              numericInput("omega_mod_weib", "Omega", .5, min = 0, max = .99, step = .1),
                              numericInput("gamma_mod_weib", "Gamma", 1, min = 0, max = 100, step = .1),
                              numericInput("kappa_mod_weib", "Kappa", 0.5, min = 0, max = 100, step = .1)
                            ),
                            fluidRow(
                              column(width = 12, class = "col-md-6",
                                     div(
                                       class = "d-flex flex-column",
                                       style = "height:90vh;",
                                       plotOutput("mod_weibullplots", height = "100%")
                                     )),
                              column(width = 12, class = "col-md-6", tags$div(style = "padding:15px;",
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
                                  tags$li("[0, 1]"),
                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                ),
                                tags$h4("Gamma"),
                                tags$ul(
                                  tags$li("[0, inf]"),
                                  tags$li("Determines the shape of the hazard rate."),
                                  tags$li("If gamma > 1, the hazard monotonically increases. If gamma < 1, hazard has a bathtub shape.")
                                ),
                                tags$h4("Kappa"),
                                tags$ul(
                                  tags$li("[0, inf]"),
                                  tags$li("Determines the shape of the hazard rate."),
                                  tags$li("If gamma < 1, the lowest point of hazard is at t_min = (sqrt(gamma) - gamma) / kappa.")
                                )
                              ))
                            )
                            )
                  ),
                  nav_panel(title="Linear-Exponential", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always", 
                              numericInput("d_lin_exp", "Duration of study", 10),
                              numericInput("f_lin_exp", "Frequency of observation", 2),
                              numericInput("omega_lin_exp", "Omega", .5, min = 0, max = .99, step = .1),
                              numericInput("gamma_lin_exp", "Gamma", 0, min = -1, max = 1, step = .1)
                            ),
                            fluidRow(
                              column(width = 12, class = "col-md-6",
                                     div(
                                       class = "d-flex flex-column",
                                       style = "height:90vh;",
                                       plotOutput("lin_expplots", height = "100%")
                                     )),
                              column(width = 12, class = "col-md-6", tags$div(style = "padding:15px;",
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
                                  tags$li("[0, 1]"),
                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                ),
                                tags$h4("Gamma"),
                                tags$ul(
                                  tags$li("[-1, 1]"),
                                  tags$li("Determines the shape of the hazard rate."),
                                  tags$li("If gamma > 0, the hazard linearly increases. If gamma < 0, the hazard linearly decreases.")
                                )
                              ))
                            )
                            )
                  ),
                  nav_panel(title="Log-Logistic", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always",  
                              numericInput("d_log", "Duration of study", 10),
                              numericInput("f_log", "Frequency of observation", 2),
                              numericInput("omega_log", "Omega", .5, min = 0, max = .99, step = .1),
                              numericInput("gamma_log", "Gamma", 1, min = 0, max = 100, step = .1)
                            ),
                            fluidRow(
                              column(width = 12, class = "col-md-6",
                                     div(
                                       class = "d-flex flex-column",
                                       style = "height:90vh;",
                                       plotOutput("logplots", height = "100%")
                                     )),
                              column(width = 12, class = "col-md-6", tags$div(style = "padding:15px;",
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
                                  tags$li("[0, 1]"),
                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                ),
                                tags$h4("Gamma"),
                                tags$ul(
                                  tags$li("[0, inf]"),
                                  tags$li("Determines the shape of the hazard rate."),
                                  tags$li("If gamma < 1, the hazard monotonically decreases. If gamma > 1, the hazard has a bell shape.")
                                )
                              ))
                            )
                            )
                  ),
                  nav_panel(title="Gompertz", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always",  
                              numericInput("d_gomp", "Duration of study", 10),
                              numericInput("f_gomp", "Frequency of observation", 2),
                              numericInput("omega_gomp", "Omega", .5, min = 0, max = .99, step = .1),
                              numericInput("gamma_gomp", "Gamma", 1, min = -100, max = 100, step = .1)
                            ),
                            fluidRow(
                              column(width = 12, class = "col-md-6",
                                      div(
                                        class = "d-flex flex-column",
                                        style = "height:90vh;",
                                        plotOutput("gompplots", height = "100%")
                                      )),
                              column(width = 12, class = "col-md-6", tags$div(style = "padding:15px;",
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
                                  tags$li("[0, 1]"),
                                  tags$li("Proportion of participants that is expected to drop out at some point during the study")
                                ),
                                tags$h4("Gamma"),
                                tags$ul(
                                  tags$li("gamma ≠ 0"),
                                  tags$li("Determines the shape of the hazard rate."),
                                  tags$li("If gamma > 0, the hazard exponentially increases. If gamma < 0, the hazard exponentially decreases.")
                                )
                              ))
                            )
                            )
                  ),
                  nav_panel(title="Non-parametric Survival function", 
                            layout_sidebar(sidebar = sidebar(
                              width = 200,   
                              open = "always",  
                              radioButtons("option", "Choose an option:", 
                                           choices = list("Enter the survival probability function, i.e., the expected proportion of participants remaining in the study for each measurement occasion" = "opt1", 
                                                          "Enter the hazard probability function, i.e., the probability of dropping out between j and j-1" = "opt2"), selected = "opt1"),
                              conditionalPanel(condition = "input.option == 'opt1'", 
                                               textInput("meas.occ", "Measurement occasions", "0, 1, 2, 5, 8, 9, 10"),
                                               textInput("survival.hand", "Expected proportion of participants remaining in the study at each measurement occasion, S(j)", "1, .95, .9, .8, .7, .65, .6")
                                               ),
                              conditionalPanel(condition = "input.option == 'opt2'",
                                               textInput("meas.occ", "Measurement occasions (comma delimited)", "0, 1, 2, 5, 8, 9, 10"),
                                               textInput("hazard.hand", "h(j), Expected risk of dropping between measurement occasion j and j-1", ".1, .15, .2, .25, .3, .5")
                                               )
                            ), 
                            fluidRow(
                              column(width = 12, class = "col-md-6",
                                      div(
                                        class = "d-flex flex-column",
                                        style = "height:90vh;",
                                        plotOutput("hand.plots", height = "100%")
                                            )),
                              column(width = 12, class = "col-md-6", 
                                     tags$h4("Choosing an option"),
                                     tags$p("Here, you can choose between specifying your own survival probability function (option 1) or your own hazard probability function (option 2). Whatever is more intuitive for you is the better option. Either way, you will see both the plot for the survival probability function as well as the plot for the hazard probability function."),
                                     tags$h4("Option 1: Specifying a survival probability function"),
                                     tags$p("First, the measurement occasions need to be specified. These can be non-equidistant, meaning that the time interval between observations does not need to be constant. If you measure, say in week 0 (baseline), in week 1, in week 5, and in week 10, then you should input '0, 1, 5, 10' for the measurement occasions. Next, you input the expected proportion of participants still remaining in the study for each of these measurement occasions (starting with 1 for the first measurement). For example, if you expect 10% of participants to drop out after the first measurement, 15% after the second measurement, and another 30% after the third measurement, you should input '1, .9, .85, .7'. Note that the first value will always be 1 as there can be no dropout before the first measurement occasion."),
                                     tags$h4("Option 2: Specifying a hazard probability function"),
                                     tags$p("The specification of the measurement occasions is the same as for option 1. Next, you input the expected risk of dropout between measurement occasion j and j - 1 for each j = 1, 2, ..., n. Note that the value for occasion j = 0 is undefined because there can be no dropout before the first measurement occasion. Using the same example measurement occasions as in option 1, if you expect the risk of dropping out in the first week to be .10, between week 1 and week 5 to be .20, and between week 5 and week 10 to be .40, your input should be '.1, .2, .4'. Note that the number of elements in expected risk need to be one less that the amount of measurement occasions as h(j=0) is undefined.")
                            )
                            )
                            )               
                            ),
                  nav_panel(title="More?", 
                            tags$h4("If you are missing information or a certain functionality in this ShinyApp, please feel free to contact me!"),
                            tags$p("Via e-mail:"),
                            tags$p("u.c.losener1@uu.nl")
                            ),
                  nav_spacer(),
                  nav_menu(
                    title = "Links",
                    nav_item(tags$a("GitHub", href = "https://github.com/ulrichlosener/BayesianSSD_attrition")),
                    nav_item(tags$a("Contact", href = "https://www.uu.nl/staff/UCLosener2/Contact"))
                  )
)

