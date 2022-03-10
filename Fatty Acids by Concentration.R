# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)
library(tidyr)

# Load data ----

FAs <- read.csv("FAs_concentration.csv", header = TRUE)
str(FAs)
FAs$ID <- as.factor(FAs$ID)
FAs$corral <- as.factor(FAs$corral)
FAs$sample.type <- as.factor(FAs$sample.type)
FAs$date <- as.Date(FAs$date,
                    format = "%b. %d, %Y")

treatments <- data.frame(corral = 
                           as.factor(c("B", "C", "D", "E", "F", "G", "H", "I")),
                         MPconcentration = as.numeric(c(0, 414, 29240, 100, 6, 
                                                        7071, 0, 1710)))

FAs <- left_join(FAs, 
                 treatments,
                 by = "corral")

# Assign corral A MP concentration of 24 for zoop_FA

FAs$MPconcentration[is.na(FAs$MPconcentration)] <- 24

# Separate out by sample type ----

perch_FA <- subset(FAs, sample.type == "perch")
zoop_FA <- subset(FAs, sample.type == "zooplankton")
food_FA <- subset(FAs, sample.type == "fish food")

# Add in perch biometrics data ----

perch_biometrics <- read.csv("perch_biometrics.csv", header = TRUE)

str(perch_biometrics)

perch_biometrics$corral <- as.factor(perch_biometrics$corral)
perch_biometrics$ID <- as.factor(perch_biometrics$ID)
perch_biometrics$sex <- as.factor(perch_biometrics$sex)

# Combine with perch FA data

perch_FA2 <- left_join(perch_FA, perch_biometrics, 
                       by = c("ID", "corral"))

# Effect of MP exposure on body weight----

perch_FA2$date2 <- 
  as.numeric(scale(as.numeric(perch_FA2$date), center = TRUE))

perchweightmod1 <- glmmTMB(body.weight ~ 
                             log(MPconcentration + 1) +
                             (1 | corral),
                           data = perch_FA2)

plotResiduals(simulateResiduals(perchweightmod1))

summary(perchweightmod1)  # no effect

perchweight_predict <- predict(perchweightmod1,
                                re.form = NA, 
                                se.fit = TRUE)

perchweight_predict$upper <- with(perchweight_predict,
                                   fit + 1.98 * se.fit)

perchweight_predict$lower <- with(perchweight_predict,
                                   fit - 1.98 * se.fit)


ggplot(perch_FA2) +
  geom_point(aes(x = MPconcentration,
                 y = body.weight)) +
  geom_ribbon(aes(x = MPconcentration,
                  ymin = perchweight_predict$lower,
                  ymax = perchweight_predict$upper),
              alpha = 0.3) +
  geom_line(aes(x = MPconcentration,
                y = perchweight_predict$fit),
            size = 1) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Body Weight (g)") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

# Add population data ----

fish_pop <- read.csv("fish_pop.csv", header = TRUE)

str(fish_pop)

fish_pop$corral <- as.factor(fish_pop$corral)
fish_pop$MP.concentration <- as.numeric(fish_pop$MP.concentration)

fish_pop$perch.surv <- with(fish_pop, YP.end/YP.start)

ggplot(fish_pop) +
  geom_point(aes(x = MP.concentration,
                 y = perch.surv)) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

# Perch Analysis----

## Total FAs ----

# Convert date to a centered numeric value

perchFAmod1 <- glmmTMB(total_FAs ~ 
                         log(MPconcentration + 1) +
                         (1 | corral),
                       data = perch_FA2)

plotResiduals(simulateResiduals(perchFAmod1))

summary(perchFAmod1)  # no effect

perch_tot_FA_predict <- predict(perchFAmod1,
                                re.form = NA, 
                                se.fit = TRUE)

perch_tot_FA_predict$upper <- with(perch_tot_FA_predict,
                                   fit + 1.98 * se.fit)

perch_tot_FA_predict$lower <- with(perch_tot_FA_predict,
                                   fit - 1.98 * se.fit)


### Plot ----

png("Perch Total FA Plot.png",
    width = 12,
    height= 8, 
    units = "cm",
    res = 600)

ggplot(perch_FA) +
  geom_point(aes(x = MPconcentration,
                 y = total_FAs,
                 fill = corral),
             size = 4,
             shape = 21) +
  geom_ribbon(aes(x = MPconcentration,
                  ymin = perch_tot_FA_predict$lower,
                  ymax = perch_tot_FA_predict$upper),
              alpha = 0.3) +
  geom_line(aes(x = MPconcentration,
                y = perch_tot_FA_predict$fit),
            size = 2) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Total fatty acids (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_fill_brewer(type = "qual",
                    palette = 1,
                    name = "Corral") +
  theme1

dev.off()

## ARA, EPA, DHA, total HUFAs -----

### ARA ----

perchARAmod1 <- glmmTMB(C_20.4n.6 ~ 
                          log(MPconcentration + 1) +
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchARAmod1))

summary(perchARAmod1)  # no effect

perchARA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchARA_pred <- predict(perchARAmod1,
                         newdata = perchARA_sim,
                         se.fit = TRUE)

perchARA_sim$pred <- perchARA_pred$fit
perchARA_sim$upper <- perchARA_pred$fit + 1.96*perchARA_pred$se.fit
perchARA_sim$lower <- perchARA_pred$fit - 1.96*perchARA_pred$se.fit

### EPA ----

perchEPAmod1 <- glmmTMB(C_20.5n.3 ~ 
                          log(MPconcentration + 1) +
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchEPAmod1))

summary(perchEPAmod1)  # no effect

perchEPA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchEPA_pred <- predict(perchEPAmod1,
                         newdata = perchEPA_sim,
                         se.fit = TRUE)

perchEPA_sim$pred <- perchEPA_pred$fit
perchEPA_sim$upper <- perchEPA_pred$fit + 1.96*perchEPA_pred$se.fit
perchEPA_sim$lower <- perchEPA_pred$fit - 1.96*perchEPA_pred$se.fit

### DHA ----

perchDHAmod1 <- glmmTMB(C_22.6n.3 ~ 
                          log(MPconcentration + 1) + 
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchDHAmod1))

summary(perchDHAmod1)  # strong effect

perchDHA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchDHA_pred <- predict(perchDHAmod1,
                         newdata = perchDHA_sim,
                         se.fit = TRUE)

perchDHA_sim$pred <- perchDHA_pred$fit
perchDHA_sim$upper <- perchDHA_pred$fit + 1.96*perchDHA_pred$se.fit
perchDHA_sim$lower <- perchDHA_pred$fit - 1.96*perchDHA_pred$se.fit

### HUFA ----

perchHUFAmod1 <- glmmTMB(HUFAs ~ 
                          log(MPconcentration + 1) + 
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchHUFAmod1))

summary(perchHUFAmod1)  # slight effect

perchHUFA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchHUFA_pred <- predict(perchHUFAmod1,
                         newdata = perchHUFA_sim,
                         se.fit = TRUE)

perchHUFA_sim$pred <- perchHUFA_pred$fit
perchHUFA_sim$upper <- perchHUFA_pred$fit + 1.96*perchHUFA_pred$se.fit
perchHUFA_sim$lower <- perchHUFA_pred$fit - 1.96*perchHUFA_pred$se.fit

### Plot ----

colours <- c("DHA" = "orange",
             "EPA" = "blue",
             "ARA" = "purple",
             "Total HUFAs" = "red")

png("Perch ARA, EPA, DHA, and HUFAs Plot.png",
    width = 12,
    height= 8, 
    units = "cm",
    res = 600)

ggplot(perch_FA) +
  geom_ribbon(data = perchARA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ARA"),
              alpha = 0.1) +
  geom_line(data = perchARA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "ARA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.4n.6,
                 colour = "ARA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = perchEPA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "EPA"),
              alpha = 0.1) +
  geom_line(data = perchEPA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "EPA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.5n.3,
                 colour = "EPA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = perchDHA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "DHA"),
              alpha = 0.1) +     
  geom_line(data = perchDHA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "DHA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3,
                 colour = "DHA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = perchHUFA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "Total HUFAs"),
              alpha = 0.1) +     
  geom_line(data = perchHUFA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "Total HUFAs")) +
  geom_point(aes(x = MPconcentration,
                 y = HUFAs,
                 colour = "Total HUFAs"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = colours,
                      name = "") +
  scale_fill_manual(values = colours,
                    name = "") +
  theme1

dev.off()

## n-3/n-6, DHA/ARA, OA/PA ----

### Put data into long form ----

perch_FA2$n3.n6 <- with(perch_FA2, total_N.3_PUFA / total_N.6_PUFA)
perch_FA2$DHA.ARA <- with(perch_FA2, C_22.6n.3 / C_20.4n.6)
perch_FA2$OA.PA <- with(perch_FA2, C_18.1n.9 / C_16.1n.7)

names(perch_FA2)

perch_FA_long <- 
  perch_FA2[,c(1:6, 17, 26, 35, 45:48, 51:63)] %>%
  pivot_longer(names(perch_FA2)[c(61:63)],
               names_to = "metric",
               values_to = "value")

perch_FA_long$metric <- as.factor(perch_FA_long$metric)

### DHA/ARA ----

perchDHA.ARAmod1 <- glmmTMB(value ~ 
                          log(MPconcentration + 1) + date2 + (1 | corral),
                          data = subset(perch_FA_long, metric == "DHA.ARA"))

plotResiduals(simulateResiduals(perchDHA.ARAmod1))

summary(perchDHA.ARAmod1)  # weak effect

perchDHA.ARA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchDHA.ARA_pred <- predict(perchDHA.ARAmod1,
                         newdata = perchDHA.ARA_sim,
                         se.fit = TRUE)

perch_FA_long$pred <- rep(NA, times = nrow(perch_FA_long))
perch_FA_long$upper <- rep(NA, times = nrow(perch_FA_long))
perch_FA_long$lower <- rep(NA, times = nrow(perch_FA_long))

perch_FA_long$pred[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit

perch_FA_long$upper[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit + 1.96*perchDHA.ARA_pred$se.fit

perch_FA_long$lower[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit - 1.96*perchDHA.ARA_pred$se.fit

### Plot ----

png("Perch FA Ratios Plot.png",
    width = 19,
    height= 15, 
    units = "cm",
    res = 600)

ggplot(perch_FA_long) +
  geom_point(aes(x = MPconcentration,
                 y = value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Ratio") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~metric,
             labeller = as_labeller(c("DHA.ARA" = "DHA:ARA",
                                      "n3.n6" = "n-3:n-6",
                                      "OA.PA" = 
                                        "Oleic acid:Palmitoleic acid")),
             ncol = 2) +
  theme1

dev.off()


# Zooplankton Analysis----

# Convert date to a centered numeric value

zoop_FA$date2 <- 
  as.numeric(scale(as.numeric(zoop_FA$date), center = TRUE))

# Convert date to a factor

zoop_FA$date2 <- as.factor(zoop_FA$date)

## Total FAs----

zoopFAmod1 <- glmmTMB(total_FAs ~
                        log(MPconcentration+1)*date2 +
                        (1 | corral),
                      data = zoop_FA)

plotResiduals(simulateResiduals(zoopFAmod1))

summary(zoopFAmod1)

zoopFAmod2 <- glmmTMB(total_FAs ~
                        MPconcentration + date2 +
                        (1 | corral),
                      data = zoop_FA)

anova(zoopFAmod1, zoopFAmod2)  # weak interaction

summary(zoopFAmod2)

AICc(zoopFAmod1, zoopFAmod2)

zoop_tot_FA_predict <- predict(zoopFAmod1,
                               re.form = NA,
                               se.fit = TRUE)

zoop_tot_FA_predict$upper <- with(zoop_tot_FA_predict,
                                   fit + 1.98 * se.fit)

zoop_tot_FA_predict$lower <- with(zoop_tot_FA_predict,
                                   fit - 1.98 * se.fit)

### Plot ----

png("Zooplankton Total FA Plot.png",
    width = 19,
    height= 8, 
    units = "cm",
    res = 600)

ggplot(zoop_FA) +
  geom_point(aes(x = MPconcentration,
                 y = total_FAs,
                 fill = corral),
             size = 4,
             shape = 21) +
  geom_ribbon(aes(x = MPconcentration,
                  ymin = zoop_tot_FA_predict$lower,
                  ymax = zoop_tot_FA_predict$upper),
              alpha = 0.3) +
  geom_line(aes(x = MPconcentration,
                y = zoop_tot_FA_predict$fit),
            size = 1) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Total fatty acids (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_fill_brewer(type = "qual",
                    palette = 1,
                    name = "Corral") +
  facet_wrap(~date) +
  theme1

dev.off()



## ARA, EPA, DHA, total HUFAs ----

### ARA ----

zoopARAmod1 <- glmmTMB(C_20.4n.6 ~ 
                          log(MPconcentration + 1)*date2 +
                          (1 | corral),
                        data = zoop_FA)

plotResiduals(simulateResiduals(zoopARAmod1))

summary(zoopARAmod1)  # no effect

zoopARA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(zoop_FA)),
                           date2 = rep(0,
                                       times = nrow(zoop_FA)),
                           corral = rep(NA,
                                        times = nrow(zoop_FA))) 

zoopARA_pred <- predict(zoopARAmod1,
                         newdata = zoopARA_sim,
                         se.fit = TRUE)

zoopARA_sim$pred <- zoopARA_pred$fit
zoopARA_sim$upper <- zoopARA_pred$fit + 1.96*zoopARA_pred$se.fit
zoopARA_sim$lower <- zoopARA_pred$fit - 1.96*zoopARA_pred$se.fit

### EPA ----

zoopEPAmod1 <- glmmTMB(C_20.5n.3 ~ 
                          log(MPconcentration + 1) +
                          date2 +
                          (1 | corral),
                        data = zoop_FA)

plotResiduals(simulateResiduals(zoopEPAmod1))

summary(zoopEPAmod1)  # no effect

zoopEPA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(zoop_FA)),
                           date2 = rep(0,
                                       times = nrow(zoop_FA)),
                           corral = rep(NA,
                                        times = nrow(zoop_FA))) 

zoopEPA_pred <- predict(zoopEPAmod1,
                         newdata = zoopEPA_sim,
                         se.fit = TRUE)

zoopEPA_sim$pred <- zoopEPA_pred$fit
zoopEPA_sim$upper <- zoopEPA_pred$fit + 1.96*zoopEPA_pred$se.fit
zoopEPA_sim$lower <- zoopEPA_pred$fit - 1.96*zoopEPA_pred$se.fit

### DHA ----

zoopDHAmod1 <- glmmTMB(C_22.6n.3 ~ 
                          log(MPconcentration + 1) + 
                          date2 +
                          (1 | corral),
                        data = zoop_FA)

plotResiduals(simulateResiduals(zoopDHAmod1))

summary(zoopDHAmod1)  # strong effect

zoopDHA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(zoop_FA)),
                           date2 = rep(0,
                                       times = nrow(zoop_FA)),
                           corral = rep(NA,
                                        times = nrow(zoop_FA))) 

zoopDHA_pred <- predict(zoopDHAmod1,
                         newdata = zoopDHA_sim,
                         se.fit = TRUE)

zoopDHA_sim$pred <- zoopDHA_pred$fit
zoopDHA_sim$upper <- zoopDHA_pred$fit + 1.96*zoopDHA_pred$se.fit
zoopDHA_sim$lower <- zoopDHA_pred$fit - 1.96*zoopDHA_pred$se.fit

### HUFA ----

zoopHUFAmod1 <- glmmTMB(HUFAs ~ 
                           log(MPconcentration + 1) + 
                           date2 +
                           (1 | corral),
                         data = zoop_FA)

plotResiduals(simulateResiduals(zoopHUFAmod1))

summary(zoopHUFAmod1)  # slight effect

zoopHUFA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                  to = 29240,
                                                  length.out = nrow(zoop_FA)),
                            date2 = rep(0,
                                        times = nrow(zoop_FA)),
                            corral = rep(NA,
                                         times = nrow(zoop_FA))) 

zoopHUFA_pred <- predict(zoopHUFAmod1,
                          newdata = zoopHUFA_sim,
                          se.fit = TRUE)

zoopHUFA_sim$pred <- zoopHUFA_pred$fit
zoopHUFA_sim$upper <- zoopHUFA_pred$fit + 1.96*zoopHUFA_pred$se.fit
zoopHUFA_sim$lower <- zoopHUFA_pred$fit - 1.96*zoopHUFA_pred$se.fit

### Plot ----

colours <- c("DHA" = "orange",
             "EPA" = "blue",
             "ARA" = "purple",
             "Total HUFAs" = "red")

png("Zooplankton ARA, EPA, DHA, and HUFAs Plot.png",
    width = 12,
    height= 8, 
    units = "cm",
    res = 600)

ggplot(zoop_FA) +
  geom_ribbon(data = zoopARA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ARA"),
              alpha = 0.1) +
  geom_line(data = zoopARA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "ARA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.4n.6,
                 colour = "ARA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = zoopEPA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "EPA"),
              alpha = 0.1) +
  geom_line(data = zoopEPA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "EPA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.5n.3,
                 colour = "EPA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = zoopDHA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "DHA"),
              alpha = 0.1) +     
  geom_line(data = zoopDHA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "DHA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3,
                 colour = "DHA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = zoopHUFA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "Total HUFAs"),
              alpha = 0.1) +     
  geom_line(data = zoopHUFA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "Total HUFAs")) +
  geom_point(aes(x = MPconcentration,
                 y = HUFAs,
                 colour = "Total HUFAs"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = colours,
                      name = "") +
  scale_fill_manual(values = colours,
                    name = "") +
  theme1

dev.off()






png("Zooplankton ARA, EPA, DHA, and HUFAs Plot.png",
    width = 12,
    height= 19, 
    units = "cm",
    res = 600)

ggplot(zoop_FA) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.4n.6,
                 colour = "ARA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_smooth(aes(x = MPconcentration,
                  y = C_20.4n.6,
                  colour = "ARA"),
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = C_20.5n.3,
                 colour = "EPA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_smooth(aes(x = MPconcentration,
                  y = C_20.5n.3,
                  colour = "EPA"),
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3,
                 colour = "DHA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_smooth(aes(x = MPconcentration,
                  y = C_22.6n.3,
                  colour = "DHA"),
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = HUFAs,
                 colour = "Total HUFAs"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_smooth(aes(x = MPconcentration,
                  y = HUFAs,
                  colour = "Total HUFAs"),
              method = "lm") +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = colours,
                      name = "") +
  facet_wrap(~date, ncol = 1) +
  theme1

dev.off()





# Fish food ----

png("Fish Food ARA, EPA and DHA Plot.png",
    width = 9,
    height= 9, 
    units = "cm",
    res = 600)

ggplot(food_FA) +
  geom_violin(aes(x = 1,
                   y = C_20.4n.6,
                   fill = "ARA"),
               size = 0.25) +
  geom_violin(aes(x = 2,
                   y = C_20.5n.3,
                   fill = "EPA"),
               size = 0.5) +
  geom_violin(aes(x = 3,
                   y = C_22.6n.3,
                   fill = "DHA"),
               size = 0.25) +
  geom_violin(aes(x = 4,
                   y = HUFAs,
                   fill = "Total HUFAs"),
               size = 0.25) +
  labs(x = "",
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_fill_manual(values = colours,
                    name = "") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

dev.off()
