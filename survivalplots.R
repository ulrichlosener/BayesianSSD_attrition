### Survival plots
library(ggplot2)
library(patchwork)


# Exponential ------------------------------------------------------------------
omega <- c(.3, .5, .8)
f <- 1
D <- 5
t.points <- seq(0, D, by=1/f)
time <- seq(0, 1, length=length(t.points))
n <- length(t.points)


remain.exp <- list()
for(i in 1:3){
  remain.exp[[i]] <- (1-omega[i])^time
}

hazard.exp <- list(NA, NA, NA)
for(j in 1:3){
  for(i in 1:(n-1)) {
      hazard.exp[[j]][i] <- round((remain.exp[[j]][i]-remain.exp[[j]][i+1])/remain.exp[[j]][i], digits=7)
  }
  hazard.exp[[j]] <- c(NA, hazard.exp[[j]])
}

surv.exp <- cbind(c(unlist(remain.exp[[1]]), unlist(remain.exp[[2]]), unlist(remain.exp[[3]])), rep(c(.3, .5, .8), each=n), rep(1:n, 3))
colnames(surv.exp) <- c("remain.exp", "omega", "occ")

p1 <- ggplot(surv.exp, aes(x=occ, y=remain.exp, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Exponential Survival Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion")

haz.exp <- cbind(c(unlist(hazard.exp[[1]]), unlist(hazard.exp[[2]]), unlist(hazard.exp[[3]])), rep(c(.3, .5, .8), each=n), rep(1:n, 3))
colnames(haz.exp) <- c("hazard", "omega", "occ")

p2 <- ggplot(haz.exp, aes(x=occ, y=hazard, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Exponential Hazard Function") +
  theme_minimal() +
  ylab("Hazard") +
  xlab("Measurement Occasion")


p1 + p2 + plot_layout(guides = "collect")

# Weibull ----------------------------------------------------------------------
omega <- c(.3, .5, .8)
gamma <- c(.5, 1, 1.5)

f <- 1
D <- 5
t.points <- seq(0, D, by=1/f)
time <- seq(0, 1, length=length(t.points))
n <- length(t.points)

dat.weib <- expand.grid(omega, gamma)
names(dat.weib) <- c("omega", "gamma")


remain.weib <- list()
for(i in 1:9){
  remain.weib[[i]] <- (1-dat.weib[i,1])^time^dat.weib[i,2]
}

hazard.weib <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA)
for(j in 1:9){
  for(i in 1:(n-1)) {
    hazard.weib[[j]][i] <- round((remain.weib[[j]][i]-remain.weib[[j]][i+1])/remain.weib[[j]][i], digits=7)
  }
  hazard.weib[[j]] <- c(NA, hazard.weib[[j]])
}

surv.weib <- cbind(unlist(remain.weib), rep(rep(c(.3, .5, .8), each=n), 3), rep(c(.5, 1, 1.5), each=n*3), rep(1:n, 9))
colnames(surv.weib) <- c("remain.weib", "omega", "gamma", "occ")

# gamma=.5
s.weib.g.5 <- ggplot(surv.weib[1:18,], aes(x=occ, y=remain.weib, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Survival Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion")

# gamma=1
s.weib.g1 <- ggplot(surv.weib[19:36,], aes(x=occ, y=remain.weib, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Survival Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion")

# gamma=1.5
s.weib.g1.5 <- ggplot(surv.weib[37:54,], aes(x=occ, y=remain.weib, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Survival Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion")

haz.weib <- cbind(unlist(hazard.weib), rep(rep(c(.3, .5, .8), each=n), 3), rep(c(.5, 1, 1.5), each=n*3), rep(1:n, 9))
colnames(haz.weib) <- c("hazard", "omega", "gamma", "occ")

# gamma=.5
h.weib.g.5 <- ggplot(haz.weib[1:18,], aes(x=occ, y=hazard, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Hazard Function") +
  theme_minimal() +
  ylab("Hazard") +
  xlab("Measurement Occasion")

# gamma=1
h.weib.g1 <- ggplot(haz.weib[19:36,], aes(x=occ, y=hazard, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Hazard Function") +
  theme_minimal() +
  ylab("Hazard") +
  xlab("Measurement Occasion")

# gamma=1.5
h.weib.g1.5 <- ggplot(haz.weib[37:54,], aes(x=occ, y=hazard, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("Hazard Function") +
  theme_minimal() +
  ylab("Hazard") +
  xlab("Measurement Occasion")


# save as PDF
pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma0.5.JPEG", width = 6, height = 4)
  s.weib.g.5 + h.weib.g.5 + plot_layout(guides = "collect")
dev.off()

pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma1.pdf", width = 6, height = 4)
  s.weib.g1 + h.weib.g1 + plot_layout(guides = "collect")
dev.off()

pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma1.5.pdf", width = 6, height = 4)
  s.weib.g1.5 + h.weib.g1.5 + plot_layout(guides = "collect")
dev.off()

# save as PNG
png(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma0.5.PNG")
  s.weib.g.5 + h.weib.g.5 + plot_layout(guides = "collect")
dev.off()

png(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma1.PNG")
  s.weib.g1 + h.weib.g1 + plot_layout(guides = "collect")
dev.off()

png(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/gamma1.5.PNG")
  s.weib.g1.5 + h.weib.g1.5 + plot_layout(guides = "collect")
dev.off()

# By hand ----------------------------------------------------------------------
t.points <- c(0,1,2,5,8,9,10)
time <- seq(0, 1, length=length(t.points))
n <- length(t.points)

remain.hand <- c(1, .8, .7, .65, .6, .58, .4)

hazard.hand <- NA
for(i in 1:(n-1)) {
  hazard.hand[i] <- round((remain.hand[i]-remain.hand[i+1])/remain.hand[i], digits=7)
}

surv.hand <- cbind(remain.hand, t.points)
colnames(surv.hand) <- c("remain", "occ")

s.hand <- ggplot(surv.hand, aes(x=occ, y=remain)) +
  geom_line() +
  ggtitle("Survival Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion") +
  scale_x_continuous(breaks = t.points)

haz.hand <- cbind(c(NA, hazard.hand), t.points)
colnames(haz.hand) <- c("hazard", "occ")

h.hand <- ggplot(haz.hand, aes(x=occ, y=hazard)) +
  geom_line() +
  ggtitle("Hazard Function") +
  theme_minimal() +
  ylab("Survival") +
  xlab("Measurement Occasion") +
  scale_x_continuous(breaks = t.points)

pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/hand.pdf")
s.hand + h.hand + plot_layout(guides = "collect")
dev.off()

png(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD_attrition/Plots/hand.PNG")
  s.hand + h.hand + plot_layout(guides = "collect")
dev.off()



