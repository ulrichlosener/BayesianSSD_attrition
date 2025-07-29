# create x variable (time)
time <- seq(0, 1, length.out = 1000)
n <- length(time)

# Create survival functions
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

# Parameters
omega <- .75
gamma_weib <- c(3, 1, .5)
gamma_mod_weib <- c(2,1,.5)
kappa_mod_weib <- .5
gamma_log <- c(2,1,.7)
gamma_lin_exp <- c(.8, 0, -.8)
gamma_gomp <- c(5, .5, -2)


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

ggplot(df_weib, aes(x = time, y = survival, color = gamma)) +
  geom_line(linewidth = 1) +
  labs(title = "Weibull Survival Function",
       x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = c(0.5, 0.15),
        legend.direction="horizontal",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"))


ggplot(df_mod_weib, aes(x = time, y = survival, color = gamma)) +
  geom_line(linewidth = 1) +
  labs(title = "Weibull Survival Function",
       x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = c(0.5, 0.15),
        legend.direction="horizontal",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"))













######################## Make hazard plots #####################################



















