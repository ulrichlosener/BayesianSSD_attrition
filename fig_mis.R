################################################################################
library(lme4)
library(data.table)
library(dplyr)

# no dropout
m <- 10000
n <- 10
BFthres <- 3
n.t.points <- 6
prop.BF1 <- rep(NA, n)
BF10 <- rep(NA, m)
a <- list()

set.seed(123)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
        a[[i]] <- getbf.mis(N=(j*10+20), t.points=c(0,1,2,3,4), dropout=F, gamma=3, omega=.2, 
                                var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=.5, 
                                fraction=1, Neff="worst", log.grow=F, hyp="h1")
        BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma = 3 --------------------------------------------------------------------

# gamma=3, omega=0.2
prop.BF1.mis1 <- rep(NA, n)

set.seed(123)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=(j*10+20), t.points=c(0,1,2,3,4), dropout=T, gamma=3, omega=.2, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=.5, 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis1[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=3, omega=0.5
prop.BF1.mis2 <- rep(NA, n)

set.seed(123)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=(j*10+20), t.points=c(0,1,2,3,4), dropout=T, gamma=3, omega=.5, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=.5, 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis2[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=3, omega=0.8
prop.BF1.mis3 <- rep(NA, n)

set.seed(123)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      tryCatch({
        a[[i]] <- getbf.mis(N=(j*10+20), t.points=c(0,1,2,3,4), dropout=T, gamma=3, omega=.8, 
                            var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=.5, 
                            fraction=1, Neff="worst", log.grow=F, hyp="h1")
        BF10[i] <- unlist(a[[i]][1])
      }, error = function(e){cat("Error :", conditionMessage(e), "\n")})
    }
  })
  prop.BF1.mis3[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# plot for gamma=3
g1 <- as.data.frame(cbind(rep(c(0, .2, .5, .8), each=n), rep(seq(30, 120, length.out=10), 4), c(prop.BF1, prop.BF1.mis1, prop.BF1.mis2, prop.BF1.mis3)))
colnames(g1) <- c("omega", "N", "power")

p1 <- ggplot(g1, aes(x=N, y=power, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("gamma=3")



pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/Plots/N.mis.pdf", width = 6, height = 4)

ggplot(g1, aes(x=N, y=power, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  theme_classic() +
  theme(legend.position = c(.8, .2)) +
  xlab("N")

dev.off()


# gamma = 1 --------------------------------------------------------------------

# gamma=1, omega=0.2
prop.BF1.mis4 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1, omega=.2, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis4[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=1, omega=0.5
prop.BF1.mis5 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1, omega=.5, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis5[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=1, omega=0.8
prop.BF1.mis6 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  c(0,1,2,3,4)  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1, omega=.8, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis6[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# plot for gamma=1
g2 <- as.data.frame(cbind(rep(c(0, .2, .5, .8), each=n), rep(1:n, 4), c(prop.BF1, prop.BF1.mis4, prop.BF1.mis5, prop.BF1.mis6)))
colnames(g2) <- c("omega", "duration", "power")

p2 <- ggplot(g2, aes(x=duration, y=power, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("gamma=1")

# gamma = 1/3 --------------------------------------------------------------------

# gamma=1/3, omega=0.2
prop.BF1.mis7 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1/3, omega=.2, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis7[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=1/3, omega=0.5
prop.BF1.mis8 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1/3, omega=.5, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis8[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# gamma=1/3, omega=0.8
prop.BF1.mis9 <- rep(NA, n)

for(j in 1:n){
  start <- Sys.time()
  suppressMessages({
    for(i in 1:m){
      a[[i]] <- getbf.mis(N=100, t.points=c(0,1,2,3,4), dropout=T, gamma=1/3, omega=.8, 
                          var.u0=.03, var.u1=.1, cov=0, var.e=.02, eff.size=(j*.1), 
                          fraction=1, Neff="worst", log.grow=F, hyp="h1")
      BF10[i] <- unlist(a[[i]][1])
    }
  })
  prop.BF1.mis9[j] <- length(BF10[BF10>BFthres])/m
  print(c(j, Sys.time()-start))
}

# plot for gamma=1/3
g3 <- as.data.frame(cbind(rep(c(0, .2, .5, .8), each=n), rep(1:n, 4), c(prop.BF1, prop.BF1.mis7, prop.BF1.mis8, prop.BF1.mis9)))
colnames(g3) <- c("omega", "duration", "power")

p3 <- ggplot(g3, aes(x=duration, y=power, color=as.factor(omega))) +
  geom_line() +
  scale_color_discrete(name="omega") +
  ggtitle("gamma=1/3")

# all together -----------------------------------------------------------------
library(gridExtra)
grid.arrange(p1, p2, p3, ncol=1)




