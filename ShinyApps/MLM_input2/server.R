################################################################################
## SHINY APP for visualizing MLM parameters - server ###########################
################################################################################

server <- function(input, output) {
  
  # contact info button
  observeEvent(input$contact_btn, {
    showModal(modalDialog(
      title = "Feel free to contact me with comments, questions, etc.!",
      HTML("<p><b>Email:</b> u.c.losener1@uu.nl</p>
            <p><b>Website:</b> <a href='https://www.uu.nl/staff/UCLosener2' target='_blank'>Utrecht University</a></p>"),
      easyClose = TRUE, # allows closing by clicking outside
      footer = modalButton("Close")
    ))
  })
  
  # number of effect sizes depends on number of conditions
  output$eff_size_inputs <- renderUI({
    req(input$n_cond)
    n_eff.size <- input$n_cond
    
    tagList(
      lapply(seq_len(n_eff.size), function(i) {
        numericInput(
          inputId = paste0("eff_size_", i),
          label = paste("Effect size for Condition", i-1),
          value = 0.4 * (i-1),
          step = 0.1
        )
      })
    )
  })
  
  # create data and plot curves
  output$responsecurves <- renderPlot({
    
    validate(
      need(input$var.e > 0, 'The residual variance needs to be larger than zero.'),
      need(input$var.u0 > 0, 'The variance in baseline scores needs to be larger than zero.'),
      need(input$var.u1 > 0, 'The variance in growth rate needs to be larger than zero.'),
      need(min(as.numeric(unlist(strsplit(input$timepoints,",")))) == 0, 'The first measurement must be taken at timepoint zero and all timepoints need to be larger than zero')
    )
    
    var.e=input$var.e
    var.u0=input$var.u0
    var.u1=input$var.u1
    covar.u01=input$covar.u01
    
    # Time points
    t <- as.numeric(unlist(strsplit(input$timepoints, ",")))
    if (input$loglinear==TRUE) {
      time <- log(t + 1)
    } else {
      time <- t
    }
    nr.points <- length(time)
    
    
    
    K <- input$n_cond  # number of conditions
    
    # Collect effect sizes (length K-1)
    eff_sizes <- sapply(seq_len(K), function(i) {
      req(input[[paste0("eff_size_", i)]])
      as.numeric(input[[paste0("eff_size_", i)]])
    })
    
    # Convert effect sizes into slope differences
    beta3_vec <- eff_sizes * sqrt(var.u1)
    
    # Baseline parameters
    int   <- input$int
    beta1 <- input$beta1

    # Create group indicator
    id    <- rep(seq_len(K), each = nr.points)
    treat <- rep(seq_len(K) - 1, each = nr.points)  # 0 = reference
    
    time_long <- rep(time, K)
    
    # Initialize outcome
    y <- numeric(length(id))
    
    for (k in seq_len(K)) {
      idx <- id == k
      y[idx] <- int + beta1*time + beta3_vec[k]*time 
    }
    
    # Random effects 
    u0 <- 0
    u1 <- 0
    resid <- 0
    
    y <- y + u0 + u1*time_long + resid
    
    dat <- data.frame(
      id    = factor(id),
      treat = treat,
      y     = y,
      time  = time_long,
      t     = rep(t, K)
    )
    
    # calculate prediction bands
    var.resp <- var.u0 + var.u1*time^2 + 2*time*covar.u01 + var.e
    
    alpha <- 1 - input$pred_width
    z <- qnorm(1 - alpha/2)
    
    dat$lower <- dat$y - z*sqrt(var.resp)
    dat$upper <- dat$y + z*sqrt(var.resp)
    
    # select conditions for which to display bands
    selected_ids <- as.numeric(input$show_pi)
    
    dat_ribbon <- dat[dat$id %in% selected_ids, ]
    
    ggplot(dat, aes(x = t, y = y, group = id)) +
      geom_line(aes(color=id), linewidth=1.5) +
      geom_point(aes(shape=id, color=id), size=4) +
      geom_ribbon(data=dat_ribbon, aes(ymin=lower, ymax=upper, group=id, fill=id),
        linetype = 0,
        alpha = 0.3
      ) +
      xlab("Time") +
      ylab("Response") +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14),        # axis tick labels
            axis.title = element_text(size = 16, face = "bold"), # axis titles
            legend.text = element_text(size = 16),      # legend labels
            ) +
      scale_color_discrete(labels = paste("Condition", 0:K)) +
      scale_shape_discrete(labels = paste("Condition", 0:K)) +
      scale_fill_discrete(labels = paste("Condition", 0:K)) +
      labs(color = NULL, shape = NULL, fill = NULL)
  })
  
  # Reactive plot title
  output$response_box <- renderUI({
    
    req(input$pred_width)
    
    perc <- round(input$pred_width * 100)
    
    box(
      title = paste0("Mean response curves and ", perc, "% prediction bands"),
      width = NULL,
      solidHeader = FALSE,
      status = "primary",
      plotOutput("responsecurves", height = "600px")
    )
  })
  
  # Selection of bands for conditions
  output$condition_selector <- renderUI({
    req(input$n_cond)
    K <- input$n_cond
    
    labels <- paste("Condition", 0:(K-1))
    
    checkboxGroupInput(
      inputId = "show_pi",
      label = HTML("<b>Show prediction band for</b>"),
      choices = setNames(seq_len(K), labels),
      selected = seq_len(K)  # default: show all
    )
  })
  
}
