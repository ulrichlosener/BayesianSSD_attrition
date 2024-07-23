library(shiny)
library(shinydashboard)
library(shinythemes)
library(plotly)
library(DT)
library(bslib)
library(ggplot2)

# ui <- dashboardPage(skin="yellow",   
# dashboardHeader(title = "Attrition Patterns",disable = F),
# dashboardSidebar(disable = T),
# dashboardBody(
#   titlePanel("Attrition Patterns in a Longitudinal Study as Modelled by a Survival Function"),
#   fluidRow(
#       column(width=6,  height = 450,
#              box(title = "Parameters", width = NULL, 
#                  numericInput("D", "Duration of study", 5),
#                  numericInput("f", "Frequency of observation", 1),
#                  numericInput("omega", "Omega (proportion of subjects who drop out at some point)", .3, min = 0, max = 1),
#                  numericInput("gamma", "Gamma (hazard rate)", 1, min = 0, max = 100),
#                  ),
#              box(title = "Create graphs",  width=NULL, background = "yellow",
#                  submitButton("Submit")
#                  )
#     
#             ),
#       column(width=12,height = 450,
#              box(title = "Attrition",width=4,height=540,
#                   plotOutput(outputId = "survivalplots", height=480))
#             )
#           )
#       )
# )

ui <- page_navbar(title ="Survival Functions", bg = "#ffd800",
                  nav_panel(title="Weibull", 
                            layout_sidebar(sidebar = sidebar(
                              numericInput("D", "Duration of study", 5),
                              numericInput("f", "Frequency of observation", 1),
                              numericInput("omega", "Omega (proportion of subjects who drop out at some point)", .3, min = 0, max = 1),
                              numericInput("gamma", "Gamma (if >1, dropout occurs more at the end of the study. if <1, dropout occurs more at the beginning of the study.)", 3, min = 0, max = 100),
                              submitButton("Create Graphs")
                            ),
                            plotOutput(outputId = "weibullplots", width="50%"),
                            textOutput(outputId = "weibullinfo"))),
                  
                  nav_panel(title="Exponential",
                            layout_sidebar(sidebar = sidebar(
                              numericInput("D", "Duration of study", 5),
                              numericInput("f", "Frequency of observation", 1),
                              numericInput("omega", "Omega (proportion of subjects who drop out at some point)", .3, min = 0, max = 1)
                            ),
                            plotOutput(outputId = "exponentialplots", width="50%"),
                            textOutput(outputId = "exponentialinfo"))
                            ),
                  nav_panel(title="Something else", p("content of page 3")),
                  nav_spacer()
)








