library(shiny)
library(ggplot2)
library(shinyWidgets)
library(shinydashboard)

fluidPage(
          # contact button
          tags$div(
            style = "position: fixed; bottom: 20px; right: 20px;",
            actionButton(
              inputId = "contact_btn",
              label = NULL,
              icon = icon("envelope"),
              style = "background-color: black; color: white; border-radius: 50%;
              width: 50px; height: 50px; padding: 0; font-size: 20px;",
              title = "Contact"
            )
          ),
          # title
          titlePanel(
            div(
              "Input Parameters for the Multilevel Model in Longitudinal Experiments",
              style = "color: #2C3E50;
                      font-size: 38px;
                      font-weight: 700;
                      text-align: center;
                      width: 100%;
                      margin-bottom: 40px;"
            ),
            windowTitle = "Input MLM Shiny App"
          ),
          setBackgroundColor(color = "#ffd800"),
          
          fluidRow(
            # column 1: number conditions and effect sizes
            column(width=3,
                   box(
                     title = "Number of conditions", width = NULL, solidHeader = FALSE, status = "primary",
                     
                     span(
                       HTML("<b>Enter the number of experimental conditions.</b>")
                     ),
                     
                     numericInput("n_cond", "", value = 2, min = 2, max = 5, step = 1),
                   ),
                   br(),
                   
                   box(title = "Effect size(s)", width = NULL, solidHeader = FALSE, status = "primary",
                     
                       withMathJax(
                         span(HTML(
                           "<b>Effect sizes of respective treatment interventions are defined as the standardized slope difference relative to the control condition. <br> The effect size of condition k is defined as \\(\\delta_k = \\beta_k \\sqrt{\\sigma^2_{u1}}\\), where \\(\\beta_k\\) is the regression coefficient of interaction between time and condition k, and \\( \\sigma^2_{u1} \\) the slope variance. An estimate may follow from the literature, expert opinion or expectations. <br> <br> As a rule of thumb, an effect size of 0.2/0.5/0.8 can be considered a small/medium/large effect. </b>"
                         ))
                       ),
                     br(), br(),
                     uiOutput("eff_size_inputs"),
                     
                   ),
                   br(),
                   
                   box(title = "Measurement occasions", width = NULL, solidHeader = FALSE, status = "primary",
                       textInput('timepoints', 'Enter time points (in increasing order and comma separated)', "0,1,2,5,8,9,10")
                   ),
                   checkboxInput("loglinear", "Use log-linear growth?", value = FALSE),
            ),
            # column 2: measurement occasions and regression parameters
            column(width=3,
                   box(
                     title = "User-specified a priori values of regression coefficients", width = NULL, solidHeader = FALSE, status = "primary",
                     # overall intercept (beta 1)
                     span(HTML(
                       "<b>Give a priori estimate of the overall intercept</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate of the initial mean score in all conditions (beta_0).",
                           icon("info-circle")))
                     ,
                     numericInput("int", "", value=0,label=NULL,step=.1),
                     # main effect time (beta 1)
                     span(HTML(
                       "<b>Give a priori estimate of the main effect of time</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate of the main effect of time (beta_1). If positive (negative), scores in all conditions will increase (decrease) over time.",
                           icon("info-circle")))
                     ,
                     numericInput("beta1", "", value=0,label=NULL,step=.1),
                   ),
                   br(),
                   
                   box(
                     title = "User-specified a priori values of (co-)variance components", width = NULL, solidHeader = FALSE, status = "primary",
                     span(HTML(
                       "<b>Give a priori estimate of residual variance</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate for the variability in residual scores. An estimate may follow from the literature, expert opinion or expectations.",
                           icon("info-circle")))
                     ,
                     numericInput("var.e", "", value=.1,label=NULL,step=.1),
                     
                     span(HTML(
                       "<b>Give a priori estimate of variance in baseline scores</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate for the variability in scores at time point zero. An estimate may follow from the literature, expert opinion or expectations.",
                           icon("info-circle")))
                     ,
                     numericInput("var.u0", "", value=.1,label=NULL,step=.1),
                     
                     
                     span(HTML(
                       "<b>Give a priori estimate of variance in growth rate</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate for the variability in growth rates. An estimate may follow from the literature, expert opinion or expectations. Note that this value affects the magnitude of the interaction between time and condition.",
                           icon("info-circle")))
                     ,
                     numericInput("var.u1", "", value=.1,label=NULL,step=.1),
                     
                     span(HTML(
                       "<b>Give a priori estimate of covariance between baseline scores and growth rate</b>"),
                       div(style = "display:inline-block;",
                           title = "User-specified a priori estimate for the association between scores at time point zero and growth rate. An estimate may follow from the literature, expert opinion or expectations.",
                           icon("info-circle")))
                     ,
                     numericInput("covar.u01", "", value=.01,label=NULL,step=.1)
                   ),
            ),
            # column 3: plot
            column(
              width = 6,
              box(
                title = "Mean response curves and 95% prediction bands",
                width = NULL,
                solidHeader = FALSE,
                status = "primary",
                plotOutput(
                  "responsecurves",
                  height = "600px"  # default ~400px
                )
              )
            )
          )
)
