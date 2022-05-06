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
# Scale date

perch_FA2$date2 <- 
  as.numeric(scale(as.numeric(perch_FA2$date), center = TRUE))

# Effect of MP exposure on body weight----

concentrations <- 
  perch_FA2 %>% 
  group_by(corral) %>% 
  summarize(MPconcentration = unique(MPconcentration)[1])

perch_biometrics2 <- left_join(perch_biometrics, concentrations, by = "corral")

perchweightmod1 <- glmmTMB(body.weight ~ 
                             log(MPconcentration + 1) +
                             (1 | corral),
                           data = subset(perch_biometrics2,
                                         !is.na(body.weight)))

plotResiduals(simulateResiduals(perchweightmod1))

summary(perchweightmod1)  # no effect

perchweight_predict <- predict(perchweightmod1,
                                re.form = NA, 
                                se.fit = TRUE)

perchweight_predict$upper <- with(perchweight_predict,
                                   fit + 1.98 * se.fit)

perchweight_predict$lower <- with(perchweight_predict,
                                   fit - 1.98 * se.fit)


ggplot(subset(perch_biometrics2,
              !is.na(body.weight))) +
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
fish_pop$MPconcentration <- as.numeric(fish_pop$MPconcentration)

fish_pop$perch.surv <- with(fish_pop, YP.end/YP.start)

ggplot(fish_pop) +
  geom_point(aes(x = MPconcentration,
                 y = perch.surv)) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

# Combine with FA data


perch_FA2 <- left_join(perch_FA2, fish_pop, 
                       by = c("corral", "MPconcentration"))

plot(perch.surv ~ total_FAs, data = perch_FA2)
plot(total_FAs ~ YP.end, data = perch_FA2)

# Perch Analysis----

## Total FAs ----

# Convert date to a centered numeric value

perchFAmod1 <- glmmTMB(total_FAs ~ 
                         log(MPconcentration + 1) +
                         body.weight +
                         (1 | corral),
                       data = perch_FA2)

plotResiduals(simulateResiduals(perchFAmod1))

summary(perchFAmod1)  # no effect

perch_tot_FA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
                                             times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perch_tot_FA_pred <- predict(perchFAmod1,
                             newdata = perch_tot_FA_sim,
                             se.fit = TRUE)

perch_tot_FA_sim$pred <- perch_tot_FA_pred$fit
perch_tot_FA_sim$upper <- perch_tot_FA_pred$fit + 1.96*perch_tot_FA_pred$se.fit
perch_tot_FA_sim$lower <- perch_tot_FA_pred$fit - 1.96*perch_tot_FA_pred$se.fit


### Plot ----

png("Perch Total FA Plot.png",
    width = 12,
    height= 8, 
    units = "cm",
    res = 600)

ggplot() +
  geom_point(data = perch_FA2,
             aes(x = MPconcentration,
                 y = total_FAs,
                 fill = corral),
             size = 4,
             shape = 21) +
  geom_ribbon(data = perch_tot_FA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper),
              alpha = 0.3) +
  geom_line(data = perch_tot_FA_sim,
            aes(x = MPconcentration,
                y = pred),
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

## Individual FAs -----

### LA ----

perchLAmod1 <- glmmTMB(C_18.2n.6 ~ 
                          log(MPconcentration + 1) +
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchLAmod1))

summary(perchLAmod1)  # no effect

perchLA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
                                             times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchLA_pred <- predict(perchLAmod1,
                         newdata = perchLA_sim,
                         se.fit = TRUE)

perchLA_sim$pred <- perchLA_pred$fit
perchLA_sim$upper <- perchLA_pred$fit + 1.96*perchLA_pred$se.fit
perchLA_sim$lower <- perchLA_pred$fit - 1.96*perchLA_pred$se.fit

### ARA ----

perchARAmod1 <- glmmTMB(C_20.4n.6 ~ 
                          log(MPconcentration + 1) +
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchARAmod1))

summary(perchARAmod1)  # no effect

perchARA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
                                             times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchARA_pred <- predict(perchARAmod1,
                         newdata = perchARA_sim,
                         se.fit = TRUE)

perchARA_sim$pred <- perchARA_pred$fit
perchARA_sim$upper <- perchARA_pred$fit + 1.96*perchARA_pred$se.fit
perchARA_sim$lower <- perchARA_pred$fit - 1.96*perchARA_pred$se.fit

### ALA ----

perchALAmod1 <- glmmTMB(C_18.3n.3 ~ 
                          log(MPconcentration + 1) +
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchALAmod1))

summary(perchALAmod1)  # no effect

perchALA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
                                             times = nrow(perch_FA2)),
                           corral = rep(NA,
                                        times = nrow(perch_FA2))) 

perchALA_pred <- predict(perchALAmod1,
                         newdata = perchALA_sim,
                         se.fit = TRUE)

perchALA_sim$pred <- perchALA_pred$fit
perchALA_sim$upper <- perchALA_pred$fit + 1.96*perchALA_pred$se.fit
perchALA_sim$lower <- perchALA_pred$fit - 1.96*perchALA_pred$se.fit

### EPA ----

perchEPAmod1 <- glmmTMB(C_20.5n.3 ~ 
                          log(MPconcentration + 1) +
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchEPAmod1))

summary(perchEPAmod1)  # no effect

perchEPA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
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
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchDHAmod1))

summary(perchDHAmod1)  # strong effect

perchDHA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
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
                          (1 | corral),
                        data = perch_FA2)

plot(simulateResiduals(perchHUFAmod1))

summary(perchHUFAmod1)  # weak effect

perchHUFA_sim <- data.frame(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = nrow(perch_FA2)),
                           date2 = rep(0,
                                       times = nrow(perch_FA2)),
                           body.weight = rep(mean(perch_FA2$body.weight),
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

colours <- c("LA" = "red", 
             "ARA" = "orange",
             "ALA" = "yellow",
             "DHA" = "green",
             "EPA" = "blue",
             "Total HUFAs" = "purple4")

png("Perch Essential FAs Plot.png",
    width = 12,
    height= 9, 
    units = "cm",
    res = 600)

ggplot(perch_FA) +
  geom_ribbon(data = perchLA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "LA"),
              alpha = 0.2) +
  geom_line(data = perchLA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "LA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_18.2n.6,
                 colour = "LA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = perchARA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ARA"),
              alpha = 0.2) +
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
  geom_ribbon(data = perchALA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ALA"),
              alpha = 0.2) +
  geom_line(data = perchALA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "ALA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_18.3n.3,
                 colour = "ALA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = perchEPA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "EPA"),
              alpha = 0.2) +
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
              alpha = 0.2) +     
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
              alpha = 0.2) +     
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

## Ratios ----

### Put data into long form ----

perch_FA2$n3.n6 <- with(perch_FA2, total_N.3_PUFA / total_N.6_PUFA)
perch_FA2$DHA.ARA <- with(perch_FA2, C_22.6n.3 / C_20.4n.6)
perch_FA2$OA.PA <- with(perch_FA2, C_18.1n.9 / C_16.1n.7)

names(perch_FA2)

perch_FA_long <- 
  perch_FA2[,c(1:6, 17, 26, 35, 45:48, 51:68)] %>%
  pivot_longer(names(perch_FA2)[c(66:68)],
               names_to = "metric",
               values_to = "value")

perch_FA_long$metric <- as.factor(perch_FA_long$metric)

### DHA/ARA ----

perchDHA.ARAmod1 <- glmmTMB(value ~ 
                          log(MPconcentration + 1) + 
                          (1 | corral),
                          data = subset(perch_FA_long, metric == "DHA.ARA"))

plotResiduals(simulateResiduals(perchDHA.ARAmod1))

summary(perchDHA.ARAmod1)  # weak effect

perch_ratio_sim <- data.frame(MPconcentration = perch_FA2$MPconcentration,
                              body.weight = rep(mean(perch_FA2$body.weight), 
                                                times = 24),
                              date2 = rep(0,
                                       times = 24),
                              corral = rep(NA,
                                           times = 24)) 

perchDHA.ARA_pred <- predict(perchDHA.ARAmod1,
                             newdata = perch_ratio_sim,
                             se.fit = TRUE,
                             re.form = NA)

perch_FA_long$pred <- rep(NA, times = nrow(perch_FA_long))
perch_FA_long$upper <- rep(NA, times = nrow(perch_FA_long))
perch_FA_long$lower <- rep(NA, times = nrow(perch_FA_long))

perch_FA_long$pred[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit

perch_FA_long$upper[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit + 1.96*perchDHA.ARA_pred$se.fit

perch_FA_long$lower[perch_FA_long$metric == "DHA.ARA"] <- 
  perchDHA.ARA_pred$fit - 1.96*perchDHA.ARA_pred$se.fit

### n-3/n-6 ----

perchn3.n6mod1 <- glmmTMB(value ~ 
                            log(MPconcentration + 1) + 
                            scale(body.weight, center = TRUE) +
                            (1 | corral),
                            data = subset(perch_FA_long, metric == "n3.n6"))

plotResiduals(simulateResiduals(perchn3.n6mod1))

summary(perchn3.n6mod1)  # no effect

perchn3.n6_pred <- predict(perchn3.n6mod1,
                           newdata = perch_ratio_sim,
                           se.fit = TRUE,
                           re.form = NA)

perch_FA_long$pred[perch_FA_long$metric == "n3.n6"] <- 
  perchn3.n6_pred$fit

perch_FA_long$upper[perch_FA_long$metric == "n3.n6"] <- 
  perchn3.n6_pred$fit + 1.96*perchn3.n6_pred$se.fit

perch_FA_long$lower[perch_FA_long$metric == "n3.n6"] <- 
  perchn3.n6_pred$fit - 1.96*perchn3.n6_pred$se.fit

### OA/PA ----

perchOA.PAmod1 <- glmmTMB(value ~ 
                            log(MPconcentration + 1) + 
                            (1 | corral),
                          data = subset(perch_FA_long, metric == "OA.PA"))

plotResiduals(simulateResiduals(perchOA.PAmod1))

summary(perchOA.PAmod1)  # no effect

perchOA.PA_pred <- predict(perchOA.PAmod1,
                           newdata = perch_ratio_sim,
                           se.fit = TRUE,
                           re.form = NA)

perch_FA_long$pred[perch_FA_long$metric == "OA.PA"] <- 
  perchOA.PA_pred$fit

perch_FA_long$upper[perch_FA_long$metric == "OA.PA"] <- 
  perchOA.PA_pred$fit + 1.96*perchOA.PA_pred$se.fit

perch_FA_long$lower[perch_FA_long$metric == "OA.PA"] <- 
  perchOA.PA_pred$fit - 1.96*perchOA.PA_pred$se.fit

### Plot ----

png("Perch FA Ratios Plot.png",
    width = 19,
    height= 7, 
    units = "cm",
    res = 600)

ggplot(perch_FA_long) +
  geom_ribbon(aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper),
              fill = "red",
              alpha = 0.3) +
  geom_line(aes(x = MPconcentration,
                y = pred)) +
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
             ncol = 3) +
  theme1

dev.off()


# Zooplankton Analysis----

# Convert date to a centered numeric value

zoop_FA$date2 <- 
  as.numeric(scale(as.numeric(zoop_FA$date), 
                   center = min(as.numeric(zoop_FA$date)),
                   scale = 74))

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

anova(zoopFAmod1, zoopFAmod2)  # not significant

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



## Individual FAs ----

### LA ----

zoopLAmod1 <- glmmTMB(C_18.2n.6 ~ 
                         log(MPconcentration + 1) * date2 +
                         (1 | corral),
                       data = zoop_FA)

plot(simulateResiduals(zoopLAmod1))

summary(zoopLAmod1)  # no effect

zoopLA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = 1000),
                           date2 = c("2021-05-27",
                                     "2021-07-06",
                                     "2021-08-09"),
                           corral = NA)

zoopLA_pred <- predict(zoopLAmod1,
                        newdata = zoopLA_sim,
                        se.fit = TRUE)

zoopLA_sim$pred <- zoopLA_pred$fit
zoopLA_sim$upper <- zoopLA_pred$fit + 1.96*zoopLA_pred$se.fit
zoopLA_sim$lower <- zoopLA_pred$fit - 1.96*zoopLA_pred$se.fit

### ARA ----

zoopARAmod1 <- glmmTMB(C_20.4n.6 ~ 
                          log(MPconcentration + 1) * date2 +
                          (1 | corral),
                        data = zoop_FA)

plot(simulateResiduals(zoopARAmod1))

summary(zoopARAmod1)  # no effect

zoopARA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = 1000),
                           date2 = c("2021-05-27",
                                     "2021-07-06",
                                     "2021-08-09"),
                           corral = NA)

zoopARA_pred <- predict(zoopARAmod1,
                         newdata = zoopARA_sim,
                         se.fit = TRUE)

zoopARA_sim$pred <- zoopARA_pred$fit
zoopARA_sim$upper <- zoopARA_pred$fit + 1.96*zoopARA_pred$se.fit
zoopARA_sim$lower <- zoopARA_pred$fit - 1.96*zoopARA_pred$se.fit

### ALA ----

zoopALAmod1 <- glmmTMB(C_18.3n.3 ~ 
                         log(MPconcentration + 1) +
                         date2 +
                         (1 | corral),
                       data = zoop_FA)

plot(simulateResiduals(zoopALAmod1))

summary(zoopALAmod1)  # no effect

zoopALA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = 1000),
                           date2 = c("2021-05-27",
                                     "2021-07-06",
                                     "2021-08-09"),
                           corral = NA)

zoopALA_pred <- predict(zoopALAmod1,
                        newdata = zoopALA_sim,
                        se.fit = TRUE)

zoopALA_sim$pred <- zoopALA_pred$fit
zoopALA_sim$upper <- zoopALA_pred$fit + 1.96*zoopALA_pred$se.fit
zoopALA_sim$lower <- zoopALA_pred$fit - 1.96*zoopALA_pred$se.fit

### EPA ----

zoopEPAmod1 <- glmmTMB(C_20.5n.3 ~ 
                          log(MPconcentration + 1) +
                          date2 +
                          (1 | corral),
                        data = zoop_FA)

plot(simulateResiduals(zoopEPAmod1))

summary(zoopEPAmod1)  # no effect

zoopEPA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = 1000),
                           date2 = c("2021-05-27",
                                     "2021-07-06",
                                     "2021-08-09"),
                           corral = NA)

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

plot(simulateResiduals(zoopDHAmod1))

summary(zoopDHAmod1)  # no effect

zoopDHA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                 to = 29240,
                                                 length.out = 1000),
                           date2 = c("2021-05-27",
                                     "2021-07-06",
                                     "2021-08-09"),
                           corral = NA)

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

plot(simulateResiduals(zoopHUFAmod1))

summary(zoopHUFAmod1)  # no effect

zoopHUFA_sim <- expand.grid(MPconcentration = seq(from = 0,
                                                  to = 29240,
                                                  length.out = 1000),
                            date2 = c("2021-05-27",
                                      "2021-07-06",
                                      "2021-08-09"),
                            corral = NA)

zoopHUFA_pred <- predict(zoopHUFAmod1,
                          newdata = zoopHUFA_sim,
                          se.fit = TRUE)

zoopHUFA_sim$pred <- zoopHUFA_pred$fit
zoopHUFA_sim$upper <- zoopHUFA_pred$fit + 1.96*zoopHUFA_pred$se.fit
zoopHUFA_sim$lower <- zoopHUFA_pred$fit - 1.96*zoopHUFA_pred$se.fit

### Plot ----

png("Zooplankton Essential FAs Plot.png",
    width = 19,
    height= 10, 
    units = "cm",
    res = 600)

ggplot(zoop_FA) +
  geom_ribbon(data = zoopLA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "LA"),
              alpha = 0.2) +
  geom_line(data = zoopLA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "LA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_18.2n.6,
                 colour = "LA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = zoopARA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ARA"),
              alpha = 0.2) +
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
  geom_ribbon(data = zoopALA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "ALA"),
              alpha = 0.2) +
  geom_line(data = zoopALA_sim,
            aes(x = MPconcentration,
                y = pred,
                colour = "ALA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_18.3n.3,
                 colour = "ALA"),
             size = 1,
             shape = 21,
             alpha = 0.75) +
  geom_ribbon(data = zoopEPA_sim,
              aes(x = MPconcentration,
                  ymin = lower,
                  ymax = upper,
                  fill = "EPA"),
              alpha = 0.2) +
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
              alpha = 0.2) +     
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
              alpha = 0.2) +     
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
  facet_wrap(~ date2) +
  theme1

dev.off()



# Fish food ----

png("Fish Food Essential FAs Plot.png",
    width = 9,
    height= 9, 
    units = "cm",
    res = 600)

ggplot(food_FA) +
  geom_violin(aes(x = 1,
                  y = C_18.2n.6,
                  fill = "LA"),
              size = 0.25) +
  geom_violin(aes(x = 2,
                   y = C_20.4n.6,
                   fill = "ARA"),
               size = 0.25) +
  geom_violin(aes(x = 3,
                  y = C_18.3n.3,
                  fill = "ALA"),
              size = 0.5) +
  geom_violin(aes(x = 4,
                   y = C_20.5n.3,
                   fill = "EPA"),
               size = 0.5) +
  geom_violin(aes(x = 5,
                   y = C_22.6n.3,
                   fill = "DHA"),
               size = 0.25) +
  geom_violin(aes(x = 6,
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
