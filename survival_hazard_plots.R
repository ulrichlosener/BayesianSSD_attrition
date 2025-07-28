library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)

# Define the survival functions
weibull <- function(omega, gamma, time){ 
  (1-omega)^(time^gamma)
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

functions <- list(
  list(name = "Weibull", 
       func = weibull, 
       params = list(gamma = c(3, 1, 0.5))),
  
  list(name = "Log-Logistic", 
       func = log_logistic, 
       params = list(gamma = c(2, 1, 0.7))),
  
  list(name = "Linear-Exponential", 
       func = linear_exponential, 
       params = list(gamma = c(0.8, 0, -0.8))),
  
  list(name = "Modified-Weibull", 
       func = mod_weibull, 
       params = list(gamma = c(2, 1, 0.5), kappa = .5)),  # kappa fixed 
  
  list(name = "Gompertz", 
       func = gompertz, 
       params = list(gamma = c(5, 0.5, -2)))
)

# Parameters
omega <- 0.75
time_values <- seq(0, 1, length.out = 100)
delta_t <- time_values[2] - time_values[1]  # Small time increment for hazard calculation

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

# 2. Data Generation (Fixed)
generate_plot_data <- function(func_info, omega, time_values, delta_t) {
  func <- func_info$func
  name <- func_info$name
  params <- func_info$params
  
  param_grid <- expand.grid(params)
  
  map_dfr(1:nrow(param_grid), function(i) {
    current_params <- as.list(param_grid[i, ])
    param_combo <- paste(names(current_params), current_params, sep="=", collapse=", ")
    
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
      param_combo = param_combo,
      stringsAsFactors = FALSE
    )
  })
}



# 3. Plotting Function (Critical - keeps gamma in legend)
create_survival_plot <- function(data, model_name) {
  is_gamma_only <- all(str_count(data$param_combo, "=") == 1)
  
  # Apply spacing rules based on model
  if (model_name == "Modified-Weibull") {
    data <- data %>%
      mutate(
        param_combo = str_replace_all(param_combo, "=", " = ")  # Add spaces around =
      )
  }
  
  p <- ggplot(data, aes(x = time, y = survival, color = param_combo)) +
    geom_line(linewidth = 1) +
    labs(title = paste(model_name, "Survival Function"),
         subtitle = expression(omega*" = 0.75"),
         x = "Time", y = "Survival Probability") +
    theme_minimal() +
    ylim(0, 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = c(0.5, 0.15),
          legend.direction="horizontal",
          legend.title = element_blank(),
          legend.background = element_rect(fill = "white", color = "black"))
  
  if (is_gamma_only) {
    p <- p + scale_color_discrete(labels = ~ paste("gamma =", sub(".*=", "", .x)))
  }
  
  return(p)
}

all_data <- map_dfr(functions, generate_plot_data, omega, time_values, delta_t)

# Split data by model and plot
plot_list <- all_data %>% 
  split(.$model) %>% 
  imap(~ create_survival_plot(.x, .y))

# Arrange in grid (2 columns)
grid.arrange(grobs = plot_list, ncol = 2)





pdf("survival_grid.pdf",
    width = 18, height = 10)
grid.arrange(grobs = plots, ncol = 2)
dev.off()



















# Create hazard plots
hazard_plots <- all_data %>%
  group_by(model) %>%
  group_map(~ {
    current_model <- .y$model
     
    p <- ggplot(.x, aes(x = time, y = hazard, color = param_combo)) +
      geom_line(linewidth = 1) +
      labs(
        title = paste(current_model, "Hazard Function"),
        # subtitle = paste("ω =", omega),
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

###
# Add parameter count to data
all_data <- all_data %>%
  mutate(
    param_count = str_count(param_combo, "=")
  )

# Create plot groups
two_param_plots <- all_data %>%
  filter(param_count == 1) %>%  # Only gamma parameter
  group_by(model) %>%
  group_map(~ create_survival_plot(.x, .y$model))

three_param_plots <- all_data %>%
  filter(param_count > 1) %>%  # Modified Weibull
  group_by(model) %>%
  group_map(~ create_survival_plot(.x, .y$model))

# Arrange plots
grid.arrange(
  arrangeGrob(grobs = two_param_plots, ncol = 2),
  arrangeGrob(grobs = three_param_plots, ncol = 1),
  nrow = 2
)

grid.arrange(grobs = survival_plots, ncol = 2)


# save plot with survival functions
pdf("survival_grid.pdf",
    width = 18, height = 10
)
grid.arrange(survival_plots[[1]], 
             survival_plots[[2]], 
             survival_plots[[3]], 
             survival_plots[[4]],
             survival_plots[[5]], 
             nrow=2)
dev.off()




grid.arrange(hazard_plots[[1]], 
             hazard_plots[[2]], 
             hazard_plots[[3]], 
             hazard_plots[[4]],
             hazard_plots[[5]], 
             nrow=2)
