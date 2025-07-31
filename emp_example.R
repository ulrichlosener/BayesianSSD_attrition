# include this as markdwon for the example -> appendix for paper

 # calculate the BFs for their results



# CONDITIONS:

# No strategy 
# Threat 
# Training 
# Fatigue
# Threat × Training 
# Threat × Fatigue
# Fatigue × Training 
# Threat × Training × Fatigue 

# we assume that:
# condition=0 is the "no strategy" condition
# condition=1 is the "addressing the threat" condition
# condition=2 is the "training effective behavior" condition
# condition=3 is the "acknowledging and combating privacy fatigue" condition


# H1: condition on severity and susceptibility


# dat <- read.csv("waveall_raw.csv")
#   
# round(dat$susc_W1,4) == round(rowMeans(cbind(dat$susceptibility_1_W1, dat$susceptibility_2_W1, dat$susceptibility_3_W1)),4)
# round(dat$susc_W2,4) == round(rowMeans(cbind(dat$susceptibility_1_W2, dat$susceptibility_2_W2, dat$susceptibility_3_W2)),4)
# round(dat$susc_W3,4) == round(rowMeans(cbind(dat$susceptibility_1_W3, dat$susceptibility_2_W3, dat$susceptibility_3_W3)),4)
# 
# 
# mean(dat$susc_W1)
# mean(dat$susc_W2)
# mean(dat$susc_W3, na.rm=T)
# 
# mean(dat_raw$susc_W1)
# mean(dat$susc_W2)
# mean(dat$susc_W3, na.rm=T)
# 
# 
# sd(dat$susc_W1)
# sd(dat$susc_W2)
# sd(dat$susc_W3, na.rm=T)

tab5 <- dat_long %>%
  group_by(wave, conditie_) %>%
  summarise(m_sus = mean(susc_),
            sd_sus = sd(susc_),
            m_sev = mean(severity_),
            sd_sev = sd(severity_))



dat %>%
  summarize(
    Time1 = sum(rowSums(!is.na(select(., ends_with("_W1")))) > 0),
    Time2 = sum(rowSums(!is.na(select(., ends_with("_W2")))) > 0),
    Time3 = sum(rowSums(!is.na(select(., ends_with("_W3")))) > 0)
  )




################################################################################
library(dplyr)
library(tidyr)


dat <- read.csv("waveall_final.csv")

dat_raw_long <- dat %>%
  select(-notnew_W1, -deelsnew_W1) %>%  # Drop problematic character columns
  # Ensure there's an ID column to track participants
  mutate(id = row_number()) %>%  
  pivot_longer(
    cols = matches("_W[1-3]r?$"),  # Matches _W1, _W2, _W3, _W1r, _W2r, _W3r
    names_to = c("variable", "wave", "suffix"),  
    names_pattern = "(.*?)_W([1-3])(r)?",  # Breaks into (variable)_W(wave)(optional r)
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  mutate(
    wave = as.integer(wave),  # Convert wave to integer
    # Combine variable name and suffix (if it exists)
    variable = ifelse(!is.na(suffix), paste0(variable, "_", suffix), variable)
  ) %>%
  select(-suffix) %>%  # Remove the temporary suffix column
  pivot_wider(
    names_from = "variable",  
    values_from = "value"
  )

# fill in missing condition values for last wave
dat_raw_long <- dat_raw_long %>%
  group_by(id) %>%
  fill(conditie_, .direction = "updown") %>%
  ungroup()

# subset data: only variables that are used in models & only conditions 0-2
dat_long <- dat_raw_long[(dat_raw_long$conditie_==1 | dat_raw_long$conditie_==2 | dat_raw_long$conditie_==3), c("id", "wave", "conditie_", "susc_", "severity_")]

# condition as factor
dat_long$conditie_ <- as.factor(dat_long$conditie_)


# add rows with NA for missing timepoints
df_complete <- dat_long %>%
  # First complete all id-wave combinations
  complete(id, wave = 1:3) %>%
  # Then fill the condition variable with previous values from the same id
  group_by(id) %>%
  fill(conditie_, .direction = "downup") %>%  # Fills both directions to handle edge cases
  ungroup()

# check dropout pattern
retention_summary <- df_complete %>%
  group_by(wave) %>%
  summarize(
    n_participants = n_distinct(id),
    prop_participants = sum(!is.na(susc_))/319,
    n_with_data = sum(!is.na(susc_))
    )  # Replace var2 with any key variable

# complete data at wave 2 ???

surviv <- c(1, 1, .586)

# plots
library(ggplot2)

# plot regression lines for each condition for wave on susceptibility
ggplot(data  = dat_long,
       aes(x = wave, y = susc_, group=conditie_, color=conditie_)) +
  geom_point(size = 1.2,
             position = "jitter",
             aes(color=conditie_)) + # "jitter" adds some random noise to differentiate the points
  geom_smooth(method = lm,
              se     = FALSE) +
  theme_minimal() +
  labs(title = "susceptibility per wave", y = "susceptibility")


# plot regression lines for each condition for wave on severity
ggplot(data  = dat_long,
             aes(x = wave, y = severity_, group=conditie_, color=conditie_)) +
  geom_point(size = 1.2,
             position = "jitter",
             aes(color=conditie_)) + # "jitter" adds some random noise to differentiate the points
  geom_smooth(method = lm,
              se     = FALSE) +
  theme_minimal() +
  labs(title = "susceptibility per wave", y = "severity")




### ANALYSES
library(lme4)
library(lmerTest)

# empty models
# mlm_empty_susc <- lmer(susc_ ~ 1 + (1|id), dat_long)
# summary(mlm_empty_susc)
# 
# mlm_empty_sev <- lmer(severity_ ~ 1 + (1|id), dat_long)
# summary(mlm_empty_sev)


mlm_susc <- lmer(susc_ ~ wave + conditie_ + wave:conditie_ + (1 + wave | id), dat_long, REML = F)
summary(mlm_susc)

mlm_sev <- lmer(severity_ ~ wave + conditie_ + wave:conditie_ + (1 + wave|id), dat_long, REML = F)
summary(mlm_sev)


### HYPOTHESIS TESTS
library(bain)

# get Neff
get_neff_mis_mv(model=mlm_susc, N=319, t.points=c(1,2,3), surviv=list(c(1, 1, .586)))

# H1_susc: wave:condition1 > 0

est_threat_susc <- mlm_susc@beta[6]
names(est_threat_susc) <- "est"
se_threat_susc <- as.matrix(vcov(mlm_susc)[6,6])

bain(x = est_threat_susc, Sigma = se_threat_susc, n = nrow(dat_long), hypothesis = "est>0")

### SSD


# effect sizes susc
# eff.size = beta / sqrt(0.000677)
# beta_cond_1 (wave) = -0.047
# beta_cond_2 = -0.022
# beta_cond_3 = -0.059

eff.sizes <-  c(-1.8, -0.8455287, -2.267)

# effect sizes sev
# eff.size = beta / sqrt(0.1336)
# beta_cond_1 (wave) = 0.028
# beta_cond_2 = -0.0533
# beta_cond_3 = -0.078

eff.sizes <- c(0.076, -0.1458, 0.21)

library(future)
library(future.apply)
library(Matrix)
library(bain)
library(dplyr)

# cov=?

BayeSSD(eta=.8, attrition=F, params=c(1, .8, .7, .6, .5, .48), 
        m=1000, t.points=c(0,1,2,3,4), var.u0=0.01, 
        var.u1=.5, var.e=.01, cov=0, eff.sizes=c(0, .5, .8), 
        BFthres=10, fraction=1, log.grow=F, seed=NULL, 
        hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.01,
        N_max=1000)










