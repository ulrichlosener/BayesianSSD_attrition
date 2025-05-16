# SIMULATION POWER MISSINGNESS THREE GROUPS

# Arguments to function
t.points <- c(0,1,2,3,4)
var.u0 <- 0.01 
var.u1 <- 0.1
cov <- 0
var.e <- .02
fraction <- 1
eff.sizes <- c(.5,.5)
hypothesis <- "a<b<c"
log.grow <- F
beta1 <- 0
BFthres <- 5

################################################################################

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
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      avg_time <- elapsed/counter
      remaining <- avg_time * (total_combos - counter)
      
      # Print progress
      cat(sprintf("\rGamma %.1f | Omega %.1f | N %d | %d/%d (%.1f%%) | Elapsed: %.0fs | Remaining: %.0fs",
                  gamma, omega, N, counter, total_combos, 
                  counter/total_combos*100, elapsed, remaining))
      flush.console()
      
      # Run single simulation
      results[[gamma_name]][[omega_name]][i] <- getpower_mis_mv(
        m = 10000 ,
        params = list(omega, gamma),
        dropout = F,
        N = N,
        t.points = t.points,
        var.u0 = var.u0,
        var.u1 = var.u1,
        cov = cov,
        var.e = var.e,
        eff.sizes = eff.sizes,
        fraction = fraction,
        log.grow = log.grow,
        beta1 = beta1,
        hypothesis = hypothesis,
        BFthres = BFthres,
        seed=123
      )$power_bf
    }
  }
}

# simulate once without dropout
for(i in seq_along(seq_N)){
  omega0[i] <- getpower_mis_mv(
    m = 10000 ,
    params = list(omega, gamma),
    dropout = FALSE,
    N = seq_N[i],
    t.points = t.points,
    var.u0 = var.u0,
    var.u1 = var.u1,
    cov = cov,
    var.e = var.e,
    eff.sizes = eff.sizes,
    fraction = fraction,
    log.grow = log.grow,
    beta1 = beta1,
    hypothesis = hypothesis,
    BFthres = BFthres,
    seed=123
  )$power_bf
  print(i)
}


# Store in long data frames
library(tidyr)
dat_gamma_0.01_wide <- as.data.frame(results$gamma_0.01)
dat_gamma_0.01_wide$N <- seq_N
dat_gamma_0.01_wide$omega_0 <- omega0
dat_gamma_0.01_wide <- dat_gamma_0.01_wide[, c("omega_0", "omega_0.1", "omega_0.3", "omega_0.5", "N")]
dat_gamma_0.01 <- gather(dat_gamma_0.01_wide, omega, power, omega_0:omega_0.5, factor_key=TRUE)

dat_gamma_0.2_wide <- as.data.frame(results$gamma_0.2)
dat_gamma_0.2_wide$N <- seq_N
dat_gamma_0.2_wide$omega_0 <- omega0
dat_gamma_0.2_wide <- dat_gamma_0.2_wide[, c("omega_0", "omega_0.1", "omega_0.3", "omega_0.5", "N")]
dat_gamma_0.2 <- gather(dat_gamma_0.2_wide, omega, power, omega_0:omega_0.5, factor_key=TRUE)

dat_gamma_1_wide <- as.data.frame(results$gamma_1)
dat_gamma_1_wide$N <- seq_N
dat_gamma_1_wide$omega_0 <- omega0
dat_gamma_1_wide <- dat_gamma_1_wide[, c("omega_0", "omega_0.1", "omega_0.3", "omega_0.5", "N")]
dat_gamma_1 <- gather(dat_gamma_1_wide, omega, power, omega_0:omega_0.5, factor_key=TRUE)

dat_gamma_5_wide <- as.data.frame(results$gamma_5)
dat_gamma_5_wide$N <- seq_N
dat_gamma_5_wide$omega_0 <- omega0
dat_gamma_5_wide <- dat_gamma_5_wide[, c("omega_0", "omega_0.1", "omega_0.3", "omega_0.5", "N")]
dat_gamma_5 <- gather(dat_gamma_5_wide, omega, power, omega_0:omega_0.5, factor_key=TRUE)


# Plots
library(ggplot2)
p_gamma_0.01 <- 
  ggplot(data = dat_gamma_0.01, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),    # Position inside bottom-right
    legend.justification = c(1, 0),            # Anchor point
    # legend.key.size = unit(2.5, "lines"),      # Large legend keys (symbols)
    # legend.text = element_text(size = 12),     # Large legend text
    # legend.title = element_text(size = 14),    # Large legend title (if exists)
    # legend.spacing.y = unit(0.2, "cm"),        # Spacing between legend items
    # legend.background = element_rect(
    #   fill = alpha("white", 0.8),             # Semi-transparent white background
    #   color = "black",                        # Border color
    #   linewidth = 0.5                         # Border thickness
    # )
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle("Gamma = 0.01")

p_gamma_0.2 <- 
  ggplot(data = dat_gamma_0.2, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),    # Position inside bottom-right
    legend.justification = c(1, 0),            # Anchor point
    # legend.key.size = unit(2.5, "lines"),      # Large legend keys (symbols)
    # legend.text = element_text(size = 12),     # Large legend text
    # legend.title = element_text(size = 14),    # Large legend title (if exists)
    # legend.spacing.y = unit(0.2, "cm"),        # Spacing between legend items
    # legend.background = element_rect(
    #   fill = alpha("white", 0.8),             # Semi-transparent white background
    #   color = "black",                        # Border color
    #   linewidth = 0.5                         # Border thickness
    # )
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle("Gamma = 0.2")

p_gamma_1 <- 
  ggplot(data = dat_gamma_1, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),   # Position inside bottom-right
    legend.justification = c(1, 0),            # Anchor point
    # legend.key.size = unit(2.5, "lines"),      # Large legend keys (symbols)
    # legend.text = element_text(size = 12),     # Large legend text
    # legend.title = element_text(size = 14),    # Large legend title (if exists)
    # legend.spacing.y = unit(0.2, "cm"),      # Spacing between legend items
    # legend.background = element_rect(
    #   fill = alpha("white", 0.8),             # Semi-transparent white background
    #   color = "black",                        # Border color
    #   linewidth = 0.5                         # Border thickness
    # )
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle("Gamma = 1")

p_gamma_5 <- 
  ggplot(data = dat_gamma_5, mapping = aes(x = N, y = power, color = omega, linetype = omega)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted")) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),   # Position inside bottom-right
    legend.justification = c(1, 0),            # Anchor point
    # legend.key.size = unit(2.5, "lines"),      # Large legend keys (symbols)
    # legend.text = element_text(size = 12),     # Large legend text
    # legend.title = element_text(size = 14),    # Large legend title (if exists)
    # legend.spacing.y = unit(0.2, "cm"),      # Spacing between legend items
    # legend.background = element_rect(
    #   fill = alpha("white", 0.8),             # Semi-transparent white background
    #   color = "black",                        # Border color
    #   linewidth = 0.5                         # Border thickness
    # )
  ) +
  scale_y_continuous(limits = c(.45, .95)) +
  ggtitle("Gamma = 5")


# Make gridplot
library(gridExtra)
library(grid)


pdf(file="./Plots/gridplot_gammas_all.pdf", width=15, height=8)
grid.arrange(p_gamma_0.01, p_gamma_0.2, p_gamma_1, p_gamma_5, nrow = 2,
             top = textGrob("Power levels per gamma and omega",gp=gpar(fontsize=20,font=3)))
dev.off()


