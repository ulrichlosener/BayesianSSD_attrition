################################################################################
# Code for simulation power and attrition with three conditions using the ######
# 'BayesSSD' package ###########################################################
################################################################################

# install BayesSSD package
if (!require("pak")) {install.packages("pak")}
pak::pak("ulrichlosener/BayesSSD")
library(BayesSSD)

# packages needed for plotting
library(ggplot2)
library(gridExtra)
library(grid)

# Arguments to function
t.points <- c(0,1,2,3,4)
var.u0 <- 0.01 
var.u1 <- 0.1
cov <- 0
var.e <- .02
fraction <- 1
eff.sizes <- c(0, .4, .8)
hypothesis <- "a<b<c"
log.grow <- F
BFthres <- 5

################################################################################
# SIMULATION 1 - BF

# Set up sequences for N, gamma, and omega
seq_N <- seq(50, 150, by=5)
seq_gamma <- c(.01, .2, 1, 5)
seq_omega <- c(.1, .3, .5)

# Initialize nested list structure
results <- list()
total_combos <- length(seq_gamma) * length(seq_omega) * length(seq_N)
counter <- 0
start_time <- Sys.time()

for(gamma in seq_gamma) {
  gamma_name <- paste0("gamma_", gamma)
  results[[gamma_name]] <- list()
  
  for(omega in seq_omega) {
    omega_name <- paste0("omega_", omega)
    results[[gamma_name]][[omega_name]] <- numeric(length(seq_N))
    names(results[[gamma_name]][[omega_name]]) <- paste0("N_", seq_N)
    
    for(i in seq_along(seq_N)) {
      N <- seq_N[i]
      
      # Update progress counter
      counter <- counter + 1
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      avg_time <- elapsed/counter
      remaining <- avg_time * (total_combos - counter)
      
      # Print progress
      cat(sprintf("\rGamma %.1f | Omega %.1f | N %d | %d/%d (%.1f%%) | Elapsed: %.0f minutes | Remaining: %.0f minutes",
                  gamma, omega, N, counter, total_combos, 
                  counter/total_combos*100, elapsed, remaining))
      flush.console()
      
      # Run single simulation
      results[[gamma_name]][[omega_name]][i] <- get_power(
        m = 10000,
        params = c(omega, gamma),
        attrition = "weibull",
        N = N,
        t.points = t.points,
        var.u0 = var.u0,
        var.u1 = var.u1,
        cov = cov,
        var.e = var.e,
        eff.sizes = eff.sizes,
        fraction = fraction,
        log.grow = log.grow,
        hypothesis = hypothesis,
        BFthres = BFthres,
        seed = 123
      )$power_bfc
    }
  }
}

# without attrition
omega0 <- rep(NA, length(seq_N))

for(i in seq_along(seq_N)){
  omega0[i] <- get_power(
    m = 10000,
    params = c(0, 1),
    attrition = F,
    N = seq_N[i],
    t.points = t.points,
    var.u0 = var.u0,
    var.u1 = var.u1,
    cov = cov,
    var.e = var.e,
    eff.sizes = eff.sizes,
    fraction = fraction,
    log.grow = log.grow,
    hypothesis = hypothesis,
    BFthres = BFthres,
    seed = 123
  )$power_bfc
  print(i/length(seq_N))
}

# Read saved simdata
dat1 <- readRDS("simresults_bfc.RDS")
omega0 <- dat1[[1]]

dat_gamma_0.01_wide <- as.data.frame(dat1[[2]])[,1:3]
dat_gamma_0.01_wide$N <- seq_N
dat_gamma_0.01_wide$omega_0 <- omega0
dat_gamma_0.01 <- gather(dat_gamma_0.01_wide, omega, power, c(omega_0, gamma_0.01.omega_0.1, gamma_0.01.omega_0.3, gamma_0.01.omega_0.5), factor_key=TRUE)
dat_gamma_0.01$omega <- rep(c("omega=0", "omega=0.1", "omega=0.3", "omega=0.5"), each=length(seq_N))

dat_gamma_0.2_wide <- as.data.frame(dat1[[2]])[,4:6]
dat_gamma_0.2_wide$N <- seq_N
dat_gamma_0.2_wide$omega_0 <- omega0
dat_gamma_0.2 <- gather(dat_gamma_0.2_wide, omega, power, c(omega_0, gamma_0.2.omega_0.1, gamma_0.2.omega_0.3, gamma_0.2.omega_0.5), factor_key=TRUE)
dat_gamma_0.2$omega <- rep(c("omega=0", "omega=0.1", "omega=0.3", "omega=0.5"), each=length(seq_N))

dat_gamma_1_wide <- as.data.frame(dat1[[2]])[,7:9]
dat_gamma_1_wide$N <- seq_N
dat_gamma_1_wide$omega_0 <- omega0
dat_gamma_1 <- gather(dat_gamma_1_wide, omega, power, c(omega_0, gamma_1.omega_0.1, gamma_1.omega_0.3, gamma_1.omega_0.5), factor_key=TRUE)
dat_gamma_1$omega <- rep(c("omega=0", "omega=0.1", "omega=0.3", "omega=0.5"), each=length(seq_N))

dat_gamma_5_wide <- as.data.frame(dat1[[2]])[,10:12]
dat_gamma_5_wide$N <- seq_N
dat_gamma_5_wide$omega_0 <- omega0
dat_gamma_5 <- gather(dat_gamma_5_wide, omega, power, c(omega_0, gamma_5.omega_0.1, gamma_5.omega_0.3, gamma_5.omega_0.5), factor_key=TRUE)
dat_gamma_5$omega <- rep(c("omega = 0", "omega = 0.1", "omega = 0.3", "omega = 0.5"), each=length(seq_N))

# Plots
# set theme according to AMPPS guidelines
ampps_theme <- theme_bw(
  base_family = "Arial",
  base_size = 9
) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(
      size = 12,
      hjust = 0.5
    ),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.box.background = element_blank(),
    axis.text.y = element_text(angle = 0)
  )

# make plots
p_gamma_0.01 <- 
  ggplot(data = dat_gamma_0.01, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),    
    legend.justification = c(1, 0),            
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle(expression(gamma == 0.01))

p_gamma_0.2 <- 
  ggplot(data = dat_gamma_0.2, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),    
    legend.justification = c(1, 0),            
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle(expression(gamma == 0.2))

p_gamma_1 <- 
  ggplot(data = dat_gamma_1, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),   
    legend.justification = c(1, 0),            
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle(expression(gamma == 1))

p_gamma_5 <- 
  ggplot(data = dat_gamma_5, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(
    values = c("solid", "longdash", "dotdash", "dotted"),
    labels = c(
      expression(omega == 0),
      expression(omega == 0.1),
      expression(omega == 0.3),
      expression(omega == 0.5)
    )
  ) +
  scale_color_discrete(
    labels = c(
      expression(omega == 0),
      expression(omega == 0.1),
      expression(omega == 0.3),
      expression(omega == 0.5)
    )
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),  
    legend.justification = c(1, 0),     
    legend.key.width = unit(2, "cm")   
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle(expression(gamma == 5))

# save plot according to AMPPS guidelines
cairo_pdf(
  file = "./LosenerFig2.pdf",
  width = 8,
  height = 6
)

grid.arrange(
  p_gamma_0.01,
  p_gamma_0.2,
  p_gamma_1,
  p_gamma_5,
  nrow = 2
)

dev.off()



###############################################################################
# SIMULATION 2 - PMP

# Set up sequences for N, gamma, and omega
seq_N <- seq(50, 150, by=5)
seq_gamma <- c(-5, -2, .01, 5)
seq_omega <- c(0, .1, .3, .5)

# Initialize nested list structure
pmp_results <- list()
total_combos <- length(seq_gamma) * length(seq_omega) * length(seq_N)
counter <- 0
start_time <- Sys.time()

for(gamma in seq_gamma) {
  gamma_name <- paste0("gamma_", gamma)
  pmp_results[[gamma_name]] <- list()
  
  for(omega in seq_omega) {
    omega_name <- paste0("omega_", omega)
    pmp_results[[gamma_name]][[omega_name]] <- numeric(length(seq_N))
    names(pmp_results[[gamma_name]][[omega_name]]) <- paste0("N_", seq_N)
    
    for(i in seq_along(seq_N)) {
      N <- seq_N[i]
      
      # Update progress counter
      counter <- counter + 1
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      avg_time <- elapsed/counter
      remaining <- avg_time * (total_combos - counter)
      
      # Print progress
      cat(sprintf("\rGamma %.1f | Omega %.1f | N %d | %d/%d (%.1f%%) | Elapsed: %.0fs | Remaining: %.0fs",
                  gamma, omega, N, counter, total_combos, 
                  counter/total_combos*100, elapsed, remaining))
      flush.console()
      
      # Run single simulation
      pmp_results[[gamma_name]][[omega_name]][i] <- get_power(
        m = 10000 ,
        params = c(omega, gamma),
        attrition = "gompertz",
        N = N,
        t.points = t.points,
        var.u0 = var.u0,
        var.u1 = var.u1,
        cov = cov,
        var.e = var.e,
        eff.sizes = eff.sizes,
        fraction = fraction,
        log.grow = log.grow,
        hypothesis = hypothesis,
        BFthres = BFthres,
        PMPthres = .9,
        seed=123
      )$power_pmp
    }
  }
}

# Read saved simdata
dat2 <- readRDS("simresults_pmp.RDS")
omega0 <- dat2[[1]]

pmp_dat_gamma_min5_wide <- as.data.frame(dat2[[2]])[,1:3]
pmp_dat_gamma_min5_wide$N <- seq_N
pmp_dat_gamma_min5_wide$omega_0 <- omega0
pmp_dat_gamma_min5 <- gather(pmp_dat_gamma_min5_wide, omega, power, c(omega_0, gamma_.5.omega_0.1, gamma_.5.omega_0.3, gamma_.5.omega_0.5), factor_key=TRUE)
pmp_dat_gamma_min5$omega <- rep(c("omega = 0", "omega = 0.1", "omega = 0.3", "omega = 0.5"), each=length(seq_N))

pmp_dat_gamma_min2_wide <- as.data.frame(dat2[[2]])[,4:6]
pmp_dat_gamma_min2_wide$N <- seq_N
pmp_dat_gamma_min2_wide$omega_0 <- omega0
pmp_dat_gamma_min2 <- gather(pmp_dat_gamma_min2_wide, omega, power, c(omega_0, gamma_.2.omega_0.1, gamma_.2.omega_0.3, gamma_.2.omega_0.5), factor_key=TRUE)
pmp_dat_gamma_min2$omega <- rep(c("omega = 0", "omega = 0.1", "omega = 0.3", "omega = 0.5"), each=length(seq_N))

pmp_dat_gamma_0.01_wide <- as.data.frame(dat2[[2]])[,7:9]
pmp_dat_gamma_0.01_wide$N <- seq_N
pmp_dat_gamma_0.01_wide$omega_0 <- omega0
pmp_dat_gamma_0.01 <- gather(pmp_dat_gamma_0.01_wide, omega, power, c(omega_0, gamma_0.01.omega_0.1, gamma_0.01.omega_0.3, gamma_0.01.omega_0.5), factor_key=TRUE)
pmp_dat_gamma_0.01$omega <- rep(c("omega = 0", "omega = 0.1", "omega = 0.3", "omega = 0.5"), each=length(seq_N))

pmp_dat_gamma_5_wide <- as.data.frame(dat2[[2]])[,10:12]
pmp_dat_gamma_5_wide$N <- seq_N
pmp_dat_gamma_5_wide$omega_0 <- omega0
pmp_dat_gamma_5 <- gather(pmp_dat_gamma_5_wide, omega, power, c(omega_0, gamma_5.omega_0.1, gamma_5.omega_0.3, gamma_5.omega_0.5), factor_key=TRUE)
pmp_dat_gamma_5$omega <- rep(c("omega = 0", "omega = 0.1", "omega = 0.3", "omega = 0.5"), each=length(seq_N))

# make plots
pmp_p_gamma_min5 <- 
  ggplot(data = pmp_dat_gamma_min5, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),    
    legend.justification = c(1, 0),            
  ) +
  scale_y_continuous(limits = c(0.8, 1)) +
  ggtitle(expression(gamma == -5))

pmp_p_gamma_min2 <- 
  ggplot(data = pmp_dat_gamma_min2, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),    
    legend.justification = c(1, 0),            
  ) +
  scale_y_continuous(limits = c(0.8, 1)) +
  ggtitle(expression(gamma == -2))

pmp_p_gamma_0.01 <- 
  ggplot(data = pmp_dat_gamma_0.01, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "none",
    legend.position.inside = c(0.95, 0.05),   
    legend.justification = c(1, 0),           
  ) +
  scale_y_continuous(limits = c(0.8, 1)) +
  ggtitle(expression(gamma == 0.01))

pmp_p_gamma_5 <- 
  ggplot(data = pmp_dat_gamma_5, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  ampps_theme +
  scale_linetype_manual(
    values = c("solid", "longdash", "dotdash", "dotted"),
    labels = c(
      expression(omega == 0),
      expression(omega == 0.1),
      expression(omega == 0.3),
      expression(omega == 0.5)
    )
  ) +
  scale_color_discrete(
    labels = c(
      expression(omega == 0),
      expression(omega == 0.1),
      expression(omega == 0.3),
      expression(omega == 0.5)
    )
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),   
    legend.justification = c(1, 0),   
    legend.key.width = unit(2, "cm")   
  ) +
  scale_y_continuous(limits = c(0.8, 1)) +
  ggtitle(expression(gamma == 5))

# save plot according to AMPPS guidelines
cairo_pdf(
  file = "./LosenerFig3.pdf",
  width = 8,
  height = 6
)

grid.arrange(
  pmp_p_gamma_min5,
  pmp_p_gamma_min2,
  pmp_p_gamma_0.01,
  pmp_p_gamma_5,
  nrow = 2
)

dev.off()

