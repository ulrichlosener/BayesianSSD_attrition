library(ggplot2)
library(gridExtra)

#################### Figure 1 & Figure 4 #######################################

# create x variable (time)
time <- seq(0, 1, length.out = 1000)
n <- length(time)

# Create survival functions
weibull <- function(omega, gamma, time){ 
  (1-omega)^(time^gamma)
}
mod_weibull <- function(omega, gamma, kappa, time){
  exp(time^gamma*exp(kappa*(time-1))*log(1-omega))
}
log_logistic <- function(omega, gamma, time){
  (1-omega)/((1-omega)+omega*time^gamma)
}
linear_exponential <- function(omega, gamma, time){
  exp((0.5*gamma+log(1-omega))*time-0.5*gamma*time^2)
}
gompertz <- function(omega, gamma, time){
  exp((log(1-omega)/(exp(gamma)-1))*(exp(gamma*time)-1))
}

# Illustrative parameters
omega <- .75
gamma_weib <- c(3, 1, .5)
gamma_mod_weib <- c(2,1.5,.5)
kappa_mod_weib <- 2
kappa_mod_weib2 <- .1
gamma_log <- c(2,1,.7)
gamma_lin_exp <- c(.9, 0, -.9)
gamma_gomp <- c(5, .5, -5)


# Create survival data for each distribution and convert to dataframes for plotting
survival_weib <- list()
for(i in 1:length(gamma_weib)){
  survival_weib[[i]] <- weibull(omega, gamma_weib[i], time)
}
df_weib <- data.frame(
  time = rep(time, length(gamma_weib)),
  survival = unlist(survival_weib),
  gamma = as.factor(rep(gamma_weib, each = length(time)))
)

survival_mod_weib <- list()
for(i in 1:length(gamma_mod_weib)){
  survival_mod_weib[[i]] <- mod_weibull(omega, gamma_mod_weib[i], kappa_mod_weib, time)
}
df_mod_weib <- data.frame(
  time = rep(time, length(gamma_mod_weib)),
  survival = unlist(survival_mod_weib),
  gamma = as.factor(rep(gamma_mod_weib, each = length(time)))
)

survival_mod_weib2 <- list()
for(i in 1:length(gamma_mod_weib)){
  survival_mod_weib2[[i]] <- mod_weibull(omega, gamma_mod_weib[i], kappa_mod_weib2, time)
}
df_mod_weib2 <- data.frame(
  time = rep(time, length(gamma_mod_weib)),
  survival = unlist(survival_mod_weib2),
  gamma = as.factor(rep(gamma_mod_weib, each = length(time)))
)

survival_log <- list()
for(i in 1:length(gamma_log)){
  survival_log[[i]] <- log_logistic(omega, gamma_log[i], time)
}
df_log <- data.frame(
  time = rep(time, length(gamma_log)),
  survival = unlist(survival_log),
  gamma = as.factor(rep(gamma_log, each = length(time)))
)

survival_lin_exp <- list()
for(i in 1:length(gamma_lin_exp)){
  survival_lin_exp[[i]] <- linear_exponential(omega, gamma_lin_exp[i], time)
}
df_lin_exp <- data.frame(
  time = rep(time, length(gamma_lin_exp)),
  survival = unlist(survival_lin_exp),
  gamma = as.factor(rep(gamma_lin_exp, each = length(time)))
)

survival_gomp <- list()
for(i in 1:length(gamma_gomp)){
  survival_gomp[[i]] <- gompertz(omega, gamma_gomp[i], time)
}
df_gomp <- data.frame(
  time = rep(time, length(gamma_gomp)),
  survival = unlist(survival_gomp),
  gamma = as.factor(rep(gamma_gomp, each = length(time)))
)

# Hazard calculation function
calculate_hazard <- function(survival_values) {
  n <- length(survival_values)
  hazard <- numeric(n)
  hazard[1:(n-1)] <- (survival_values[1:(n-1)] - survival_values[2:n]) / survival_values[1:(n-1)]
  hazard[n] <- NA  # Last point can't be calculated
  return(hazard)
}

# Calculate hazards for each distribution and convert to dataframes for plotting
hazard_weib <- list()
for(i in 1:length(gamma_weib)){
  hazard_weib[[i]] <- calculate_hazard(survival_weib[[i]])
}
df_weib_haz <- data.frame(
  time = rep(time, length(gamma_weib)),
  survival = unlist(hazard_weib),
  gamma = as.factor(rep(gamma_weib, each = length(time)))
)

hazard_mod_weib <- list()
for(i in 1:length(gamma_mod_weib)){
  hazard_mod_weib[[i]] <- calculate_hazard(survival_mod_weib[[i]])
}
df_mod_weib_haz <- data.frame(
  time = rep(time, length(gamma_mod_weib)),
  survival = unlist(hazard_mod_weib),
  gamma = as.factor(rep(gamma_mod_weib, each = length(time)))
)

hazard_mod_weib2 <- list()
for(i in 1:length(gamma_mod_weib)){
  hazard_mod_weib2[[i]] <- calculate_hazard(survival_mod_weib2[[i]])
}
df_mod_weib_haz2 <- data.frame(
  time = rep(time, length(gamma_mod_weib)),
  survival = unlist(hazard_mod_weib2),
  gamma = as.factor(rep(gamma_mod_weib, each = length(time)))
)

hazard_log <- list()
for(i in 1:length(gamma_log)){
  hazard_log[[i]] <- calculate_hazard(survival_log[[i]])
}
df_log_haz <- data.frame(
  time = rep(time, length(gamma_log)),
  survival = unlist(hazard_log),
  gamma = as.factor(rep(gamma_log, each = length(time)))
)

hazard_lin_exp <- list()
for(i in 1:length(gamma_lin_exp)){
  hazard_lin_exp[[i]] <- calculate_hazard(survival_lin_exp[[i]])
}
df_lin_exp_haz <- data.frame(
  time = rep(time, length(gamma_lin_exp)),
  survival = unlist(hazard_lin_exp),
  gamma = as.factor(rep(gamma_lin_exp, each = length(time)))
)

hazard_gomp <- list()
for(i in 1:length(gamma_gomp)){
  hazard_gomp[[i]] <- calculate_hazard(survival_gomp[[i]])
}
df_gomp_haz <- data.frame(
  time = rep(time, length(gamma_gomp)),
  survival = unlist(hazard_gomp),
  gamma = as.factor(rep(gamma_gomp, each = length(time)))
)


######################## Make survival plots ###################################

p_weib <- ggplot(df_weib, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
            geom_line(linewidth = 1.25) +
            labs(title = "Weibull Survival Function",
                 x = "Time (proportion)", y = "Survival Probability",
                 color=expression(gamma), linetype=expression(gamma)) +
            ylim(c(0,1)) +
            theme_minimal() +
            theme(legend.position = c(0.5, 0.15),
                  legend.background = element_rect(fill = "white", color = "black"),
                  legend.direction="horizontal",
                  legend.title = element_text(size=20),
                  legend.key.size = unit(1.4,"line"),
                  title = element_text(face="bold"),
                  axis.title = element_text(face="plain")) 

p_mod_weib <- ggplot(df_mod_weib, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                geom_line(linewidth = 1.25) +
                labs(title = "Modified Weibull Survival Function",
                     subtitle = expression(kappa*" = 2"),
                     x = "Time (proportion)", y = "Survival Probability",
                     color=expression(gamma), linetype=expression(gamma)) +
                ylim(c(0,1)) +
                theme_minimal() +
                theme(legend.position = c(0.5, 0.15),
                      legend.direction="horizontal",
                      legend.title = element_text(size=20),
                      legend.background = element_rect(fill = "white", color = "black"),
                      legend.key.size = unit(1.4,"line"),
                      title = element_text(face="bold"),
                      plot.subtitle = element_text(size=15),
                      axis.title = element_text(face="plain"))
              
p_mod_weib2 <- ggplot(df_mod_weib2, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                  geom_line(linewidth = 1.25) +
                  labs(title = "Modified Weibull Survival Function",
                       subtitle = expression(kappa*" = 0.1"),
                       x = "Time (proportion)", y = "Survival Probability",
                       color=expression(gamma), linetype=expression(gamma)) +
                  ylim(c(0,1)) +
                  theme_minimal() +
                  theme(legend.position = c(0.5, 0.15),
                        legend.direction="horizontal",
                        legend.title = element_text(size=20),
                        legend.background = element_rect(fill = "white", color = "black"),
                        legend.key.size = unit(1.4,"line"),
                        title = element_text(face="bold"),
                        plot.subtitle = element_text(size=15),
                        axis.title = element_text(face="plain"))

p_log <- ggplot(df_log, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
          geom_line(linewidth = 1.25) +
          labs(title = "Log-Logistic Survival Function",
               x = "Time (proportion)", y = "Survival Probability",
               color=expression(gamma), linetype=expression(gamma)) +
          ylim(c(0,1)) +
          theme_minimal() +
          theme(legend.position = c(0.5, 0.15),
                legend.background = element_rect(fill = "white", color = "black"),
                legend.direction="horizontal",
                legend.title = element_text(size=20),
                legend.key.size = unit(1.4,"line"),
                title = element_text(face="bold"),
                axis.title = element_text(face="plain")) 
        
p_lin_exp <- ggplot(df_lin_exp, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
              geom_line(linewidth = 1.25) +
              labs(title = "Linear-Exponential Survival Function",
                   x = "Time (proportion)", y = "Survival Probability",
                   color=expression(gamma), linetype=expression(gamma)) +
              ylim(c(0,1)) +
              theme_minimal() +
              theme(legend.position = c(0.5, 0.15),
                    legend.background = element_rect(fill = "white", color = "black"),
                    legend.direction="horizontal",
                    legend.title = element_text(size=20),
                    legend.key.size = unit(1.4,"line"),
                    title = element_text(face="bold"),
                    axis.title = element_text(face="plain")) 
            
p_gomp <- ggplot(df_gomp, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
            geom_line(linewidth = 1.25) +
            labs(title = "Gompertz Survival Function",
                 x = "Time (proportion)", y = "Survival Probability",
                 color=expression(gamma), linetype=expression(gamma)) +
            ylim(c(0,1)) +
            theme_minimal() +
            theme(legend.position = c(0.5, 0.15),
                  legend.background = element_rect(fill = "white", color = "black"),
                  legend.direction="horizontal",
                  legend.title = element_text(size=20),
                  legend.key.size = unit(1.4,"line"),
                  title = element_text(face="bold"),
                  axis.title = element_text(face="plain")) 


pdf("survival_grid2.pdf", width = 8, height=12)
  grid.arrange(p_weib, p_mod_weib, p_mod_weib2, p_log, p_lin_exp, p_gomp, nrow=3)
dev.off()


######################## Make hazard plots #####################################

p_weib_haz <- ggplot(df_weib_haz, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                  geom_line(linewidth = 1.25) +
                  labs(title = "Weibull Hazard Function",
                       x = "Time (proportion)", y = "Hazard Probability",
                       color=expression(gamma), linetype=expression(gamma)) +
                  theme_minimal() +
                  theme(legend.position = c(0.5, 0.8),
                        legend.background = element_rect(fill = "white", color = "black"),
                        legend.direction="horizontal",
                        legend.title = element_text(size=20),
                        legend.key.size = unit(1.4,"line"),
                        title = element_text(face="bold"),
                        axis.title = element_text(face="plain")) +
                  ylim(c(0, .005))


p_mod_weib_haz <- ggplot(df_mod_weib_haz, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                    geom_line(linewidth = 1.25) +
                    labs(title = "Modified Weibull Hazard Function",
                         subtitle = expression(kappa*" = 2"),
                         x = "Time (proportion)", y = "Hazard Probability",
                         color=expression(gamma), linetype=expression(gamma)) +
                    theme_minimal() +
                    theme(legend.position = c(0.5, 0.8),
                          legend.direction="horizontal",
                          legend.title = element_text(size=20),
                          legend.background = element_rect(fill = "white", color = "black"),
                          legend.key.size = unit(1.4,"line"),
                          title = element_text(face="bold"),
                          plot.subtitle = element_text(size=15),
                          axis.title = element_text(face="plain"))

p_mod_weib_haz2 <- ggplot(df_mod_weib_haz2, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                    geom_line(linewidth = 1.25) +
                    labs(title = "Modified Weibull Hazard Function",
                         subtitle = expression(kappa*" = 0.1"),
                         x = "Time (proportion)", y = "Hazard Probability",
                         color=expression(gamma), linetype=expression(gamma)) +
                    theme_minimal() +
                    theme(legend.position = c(0.5, 0.8),
                          legend.direction="horizontal",
                          legend.title = element_text(size=20),
                          legend.background = element_rect(fill = "white", color = "black"),
                          legend.key.size = unit(1.4,"line"),
                          title = element_text(face="bold"),
                          plot.subtitle = element_text(size=15),
                          axis.title = element_text(face="plain")) +
                    ylim(c(0, .004))

p_log_haz <- ggplot(df_log_haz, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
              geom_line(linewidth = 1.25) +
              labs(title = "Log-Logistic Hazard Function",
                   x = "Time (proportion)", y = "Hazard Probability",
                   color=expression(gamma), linetype=expression(gamma)) +
              theme_minimal() +
              theme(legend.position = c(0.5, 0.8),
                    legend.background = element_rect(fill = "white", color = "black"),
                    legend.direction="horizontal",
                    legend.title = element_text(size=20),
                    legend.key.size = unit(1.4,"line"),
                    title = element_text(face="bold"),
                    axis.title = element_text(face="plain")) +
              ylim(c(0,.005))

p_lin_exp_haz <- ggplot(df_lin_exp_haz, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                  geom_line(linewidth = 1.25) +
                  labs(title = "Linear-Exponential Hazard Function",
                       x = "Time (proportion)", y = "Hazard Probability",
                       color=expression(gamma), linetype=expression(gamma)) +
                  theme_minimal() +
                  theme(legend.position = c(0.5, 0.09),
                        legend.background = element_rect(fill = "white", color = "black"),
                        legend.direction="horizontal",
                        legend.title = element_text(size=20),
                        legend.key.size = unit(1.4,"line"),
                        title = element_text(face="bold"),
                        axis.title = element_text(face="plain")                        ) 

p_gomp_haz <- ggplot(df_gomp_haz, aes(x = time, y = survival, color = gamma, linetype=gamma)) +
                geom_line(linewidth = 1.25) +
                labs(title = "Gompertz Hazard Function",
                     x = "Time (proportion)", y = "Hazard Probability",
                     color=expression(gamma), linetype=expression(gamma)) +
                theme_minimal() +
                theme(legend.position = c(0.5, 0.8),
                      legend.background = element_rect(fill = "white", color = "black"),
                      legend.direction="horizontal",
                      legend.title = element_text(size=20),
                      legend.key.size = unit(1.4,"line"),
                      title = element_text(face="bold"),
                      axis.title = element_text(face="plain")) +
                ylim(0,.0075)

pdf("hazard_grid2.pdf", width = 8, height=12)
  grid.arrange(p_weib_haz, p_mod_weib_haz, p_mod_weib_haz2, p_log_haz, p_lin_exp_haz, p_gomp_haz)
dev.off()
















