# Plots of different survival distributions
library(ggplot2)

x_range <- seq(0, 1, .01)
surv_weib <- weibull(omega=.5, gamma=.5, time=x_range)
dat_weib <- 

ggplot(data)









library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(gridExtra)

# Define the survival functions
weibull <- function(omega, gamma, time){ 
  (1-omega)^(time^gamma)
}
exponential <- function(omega, time){
  (1-omega)^time
}
log_logistic <- function(omega, gamma, time){
  (1-omega)/((1-omega)+omega*time^gamma)
}
linear_exponential <- function(omega, gamma, time){
  exp((0.5*gamma+log(1-omega))*time-0.5*gamma*time^2)
}
mod_weibull <- function(omega, gamma, kappa, time){
  exp(time^gamma*exp(kappa*(time-1))*log(1-omega))
}
gompertz <- function(omega, gamma, time){
  exp((log(1-omega)/(exp(gamma)-1))*(exp(gamma*time)-1))
}

# Parameters
omega <- 0.5
time_values <- seq(0, 1, length.out = 100)
delta_t <- time_values[2] - time_values[1]  # Small time increment for hazard calculation

# Create a list of functions with their custom parameter requirements
functions <- list(
  list(name = "Weibull", 
       func = weibull, 
       params = list(gamma = c(3, 1, 0.3))),
  
  list(name = "Exponential", 
       func = exponential, 
       params = list()),
  
  list(name = "Log-Logistic", 
       func = log_logistic, 
       params = list(gamma = c(2, 1, 0.2))),
  
  list(name = "Linear-Exponential", 
       func = linear_exponential, 
       params = list(gamma = c(0.5, 0, -0.8))),
  
  list(name = "Modified-Weibull", 
       func = mod_weibull, 
       params = list(gamma = c(2, 1, 0.1), kappa = 2)),  # kappa fixed 
  
  list(name = "Gompertz", 
       func = gompertz, 
       params = list(gamma = c(5, 0.5, -2)))
)

# Function to calculate hazard rate
calculate_hazard <- function(survival_values, delta_t) {
  n <- length(survival_values)
  hazard <- numeric(n)
  hazard[1:(n-1)] <- (survival_values[1:(n-1)] - survival_values[2:n]) / survival_values[1:(n-1)]
  hazard[n] <- hazard[n-1]  # Extend the last value
  hazard[is.infinite(hazard)] <- NA  # Handle division by zero
  hazard[is.nan(hazard)] <- NA       # Handle 0/0 cases
  return(hazard)
}

# Function to generate data for plotting
generate_plot_data <- function(func_info, omega, time_values, delta_t) {
  func <- func_info$func
  name <- func_info$name
  params <- func_info$params
  
  # Get all combinations of parameters
  param_grid <- expand.grid(params)
  
  # If no additional parameters (just exponential)
  if (ncol(param_grid) == 0) {
    survival <- sapply(time_values, function(t) func(omega, t))
    hazard <- calculate_hazard(survival, delta_t)
    
    return(data.frame(
      time = time_values,
      survival = survival,
      hazard = hazard,
      model = name,
      param_combo = "none"
    ))
  }
  
  # For functions with parameters
  results <- map_dfr(1:nrow(param_grid), function(i) {
    current_params <- as.list(param_grid[i, ])
    param_combo <- paste(names(current_params), current_params, sep = "=", collapse = ", ")
    
    survival <- sapply(time_values, function(t) {
      args <- c(omega = omega, current_params, time = t)
      do.call(func, args)
    })
    
    hazard <- calculate_hazard(survival, delta_t)
    
    data.frame(
      time = time_values,
      survival = survival,
      hazard = hazard,
      model = name,
      param_combo = param_combo
    )
  })
  
  return(results)
}

# Generate data for all functions
all_data <- map_dfr(functions, generate_plot_data, omega, time_values, delta_t)

# Create survival plots
survival_plots <- all_data %>%
  group_by(model) %>%
  group_map(~ {
    current_model <- .y$model
    
    p <- ggplot(.x, aes(x = time, y = survival, color = param_combo)) +
      geom_line(linewidth = 1) +
      labs(
        title = paste(current_model, "Survival Function"),
        subtitle = paste("ω =", omega),
        x = "Time",
        y = "Survival Probability",
        color = "Parameters"
      ) +
      theme_minimal() +
      ylim(0, 1) +
      scale_color_brewer(palette = "Dark2") +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12, face = "bold"))
    
    if (n_distinct(.x$param_combo) == 1) {
      p <- p + theme(legend.position = "none")
    }
    return(p)
  })

# Create hazard plots
hazard_plots <- all_data %>%
  group_by(model) %>%
  group_map(~ {
    current_model <- .y$model
    
    p <- ggplot(.x, aes(x = time, y = hazard, color = param_combo)) +
      geom_line(linewidth = 1) +
      labs(
        title = paste(current_model, "Hazard Function"),
        subtitle = paste("ω =", omega),
        x = "Time",
        y = "Hazard Rate",
        color = "Parameters"
      ) +
      theme_minimal() +
      scale_color_brewer(palette = "Dark2") +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12, face = "bold"))
    
    if (n_distinct(.x$param_combo) == 1) {
      p <- p + theme(legend.position = "none")
    }
    return(p)
  })

# Display plots (they will appear in RStudio's plot pane)
# To view specific plots, you can access them from survival_plots or hazard_plots lists
survival_plots[[6]]
hazard_plots[[5]]


grid.arrange(survival_plots[[1]], 
             survival_plots[[2]], 
             survival_plots[[3]], 
             survival_plots[[4]],
             survival_plots[[5]], 
             survival_plots[[6]],
             nrow=2)


grid.arrange(hazard_plots[[1]], 
             hazard_plots[[2]], 
             hazard_plots[[3]], 
             hazard_plots[[4]],
             hazard_plots[[5]], 
             hazard_plots[[6]],
             nrow=2)




