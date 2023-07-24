# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)
library(tidyr)
library(vegan)
library(emmeans)
library(ggeffects)
library(cowplot)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 7,
                        family = "serif"),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  )

FA.names <- c("12:0", "13:0", "14:0", "15:0", "16:0", "17:0", "18:0", "20:0", 
              "22:0", "24:0", "12:1", "14:1", "16:1n-7", "16:1n-9", "18:1n-7",
              "18:1n-9", "20:1n-9", "22:1n-9", "18:2n-6", "18:3n-6", "20:2n-6",
              "20:3n-6", "20:4n-6", "22:2n-6", "22:4n-6", "22:5n-6", "18:3n-3",
              "18:4n-3", "20:3n-3", "20:4n-3", "20:5n-3", "22:5n-3", "22:6n-3",
              "24:5n-3", "24:6n-3")

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

## Add diet info ----

diet_ca_cn$corral <- as.factor(diet_ca_cn$corral)

perch_FA3 <-
  left_join(perch_FA2,
            diet_ca_cn[,c(1:3)],
            by = "corral")

# Calculate ratios

perch_FA3$ALA.LIN <- with(perch_FA3, C_18.3n.3 / C_18.2n.6)
perch_FA3$PUFA.SMUFA <- with(perch_FA3, PUFAs / (total_SFAs + total_MUFAs))
perch_FA3$ARA.LIN <- with(perch_FA3, C_20.4n.6 / C_18.2n.6)
perch_FA3$DHA.ALA <- with(perch_FA3, C_22.6n.3 / C_18.3n.3)
perch_FA3$DHA.ARA <- with(perch_FA3, C_22.6n.3 / C_20.4n.6)

names(perch_FA3)

# Reduce dataset to fish with sex and gonad weight and diet cca data

perch_FA2sex <- 
  perch_FA3 %>% 
  filter(!is.na(sex) &
           !is.na(gonad.weight) &
           !is.na(CCA1))

nrow(perch_FA2sex)  # 18 fish remaining

## Total FAs ----

plot(total_FAs ~ TL, data = perch_FA2)

perchFAmod1 <- glmmTMB(log(total_FAs) ~ 
                         body.weight +
                         log(MPconcentration + 6) +
                         (1 | corral),
                       data = perch_FA3)

plot(simulateResiduals(perchFAmod1))

summary(perchFAmod1) 
# positive correlation with length, negative with gonad weight

anova(perchFAmod1, update(perchFAmod1, ~-corral))  # corral not significant

## Relationship of all FAs with total FAs ----

### Put data into long form

perch_FAlong <- 
  perch_FA3 %>% 
  pivot_longer(cols = c(7:16,18:25,27:34,36:44),
               names_to = "FAs")

perch_FAlong$FAs <- as.factor(perch_FAlong$FAs)

## Plot ----

# total FAs
ggplot(perch_FAlong) +
  geom_smooth(aes(x = total_FAs,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = total_FAs,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

# body size
ggplot(perch_FAlong) +
  geom_smooth(aes(x = TL,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = TL,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

# gonad size
ggplot(perch_FAlong) +
  geom_smooth(aes(x = gonad.weight,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = gonad.weight,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

`# CCA1
ggplot(perch_FAlong) +
  geom_smooth(aes(x = CCA1,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = CCA1,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

# CCA2
ggplot(perch_FAlong) +
  geom_smooth(aes(x = CCA2,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = CCA2,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

# Sample date
ggplot(perch_FAlong) +
  geom_smooth(aes(x = date,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = date,
                 y = value,
                 colour = FAs)) +
  facet_wrap(~FAs) +
  theme1

## Individual FAs -----

### Treatment effect on HUFAs ----

#### EPA ----

perchEPAmod1 <- 
  glmmTMB(C_20.5n.3 ~
            log(MPconcentration + 6) + 
            body.weight +
            (1 | corral),
          data = perch_FA3)

plot(simulateResiduals(perchEPAmod1))

summary(perchEPAmod1)  
# positively correlated with total length

perchEPA_pred <- 
  ggemmeans(perchEPAmod1,
            terms = c("body.weight")) %>%
  rename(body.weight = x)

png("Perch EPA Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA3) +
  geom_ribbon(data = perchEPA_pred,
              aes(x = body.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchEPA_pred,
            aes(x = body.weight,
                y = predicted)) +
  geom_point(aes(x = body.weight,
                 y = C_20.5n.3),
             size = 1) +
  scale_colour_viridis_d(option = "plasma",
                         name = "Total Length (cm)") +
  labs(x = "Body Weight (g)",
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  theme1

dev.off()

#### DHA ----

perchDHAmod1 <- 
  glmmTMB(C_22.6n.3 ~
            log(MPconcentration + 6) + 
            body.weight +
            (1 | corral),
          data = perch_FA3)

plot(simulateResiduals(perchDHAmod1))

summary(perchDHAmod1)  
# positively correlated with MP concentration

perchDHA_pred <- 
  ggemmeans(perchDHAmod1,
            terms = c("MPconcentration [0.001:29240, by = 200]")) %>%
  rename(MPconcentration = x)

png("Perch DHA Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA3) +
  geom_ribbon(data = perchDHA_pred,
              aes(x = MPconcentration,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchDHA_pred,
            aes(x = MPconcentration,
                y = predicted)) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3),
             size = 1) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0,1,10,100,1000,10000)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  theme1

dev.off()

#### ARA ----

perchARAmod1 <- 
  glmmTMB(C_20.4n.6 ~
            body.weight + 
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA3)

plot(simulateResiduals(perchARAmod1))

summary(perchARAmod1) 

### Effects of biometrics ----

#### ALA ----

perchALAmod1 <- 
  glmmTMB(C_18.3n.3 ~
            log(gonad.weight) + 
            body.weight +
            sex +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchALAmod1))

summary(perchALAmod1)  
# positively correlated body weight and negatively correlated with gonad weight

perchALA_pred <- 
  ggemmeans(perchALAmod1,
            terms = c("gonad.weight",
                      "body.weight")) %>%
  rename(gonad.weight = x,
         body.weight = group)

ALAplot <-
  ggplot(perch_FA2sex) +
    geom_ribbon(data = perchALA_pred,
                aes(x = gonad.weight,
                    ymin = conf.low,
                    ymax = conf.high,
                    fill = body.weight),
                alpha = 0.3) +
    geom_line(data = perchALA_pred,
              aes(x = gonad.weight,
                  y = predicted,
                  colour = body.weight)) +
    geom_point(aes(x = gonad.weight,
                   y = C_18.3n.3),
               size = 1) +
    scale_colour_viridis_d(option = "plasma",
                           name = "Body Weight (g)") +
    scale_fill_viridis_d(option = "plasma",
                         "Body Weight (g)") +
    labs(x = "Gonad Weight (g)",
         y = expression(paste("Concentration (mg "~g^-1*")"))) +
    theme1

#### LIN ----

perchLINmod1 <- 
  glmmTMB(C_18.2n.6 ~
            log(gonad.weight) + 
            body.weight +
            sex +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchLINmod1))

summary(perchLINmod1)  
# weak negative correlation with gonad weight

perchLIN_pred <- 
  ggemmeans(perchLINmod1,
            terms = c("gonad.weight")) %>%
  rename(gonad.weight = x)

LINplot <-
  ggplot(perch_FA2sex) +
    geom_ribbon(data = perchLIN_pred,
                aes(x = gonad.weight,
                    ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.3) +
    geom_line(data = perchLIN_pred,
              aes(x = gonad.weight,
                  y = predicted)) +
    geom_point(aes(x = gonad.weight,
                   y = C_18.2n.6),
               size = 1) +
    scale_colour_viridis_d(option = "plasma",
                           name = "Total Length (cm)") +
    scale_fill_viridis_d(option = "plasma",
                         "Total Length (cm)") +
    labs(x = "Gonad Weight (g)",
         y = expression(paste("Concentration (mg "~g^-1*")"))) +
    theme1

#### ALA and LIN plot ----

png("Perch ALA and LIN Plot.png",
    width = 8.84,
    height= 10, 
    units = "cm",
    res = 500)

plot_grid(ALAplot, LINplot, ncol = 1,
          align = "v",
          axis = "r",
          labels = c("A", "B"))

dev.off()

#### DHA ----

perchDHAmod2 <- 
  glmmTMB(C_22.6n.3 ~
            body.weight +
            gonad.weight +
            sex +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchDHAmod1))

summary(perchDHAmod2)  
# weak positive correlation with body weight

perchDHA_pred2 <- 
  ggemmeans(perchDHAmod2,
            terms = c("body.weight")) %>%
  rename(body.weight = x)

png("Perch DHA Plot 2.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchDHA_pred2,
              aes(x = body.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchDHA_pred2,
            aes(x = body.weight,
                y = predicted)) +
  geom_point(aes(x = body.weight,
                 y = C_22.6n.3),
             size = 1) +
  labs(x = "Body Weight (g)",
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  theme1

dev.off()

#### n-6 DPA ----

perchn6DPAmod1 <- 
  glmmTMB(log(C_22.5n.6) ~
            gonad.weight + 
            body.weight +
            sex +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchn6DPAmod1))

summary(perchn6DPAmod1)  
# no effects

## Ratios ----

### ALA/LIN (diet) ----

plot(ALA.LIN ~ corral, data = perch_FA2sex)

perchALA.LINmod1 <- 
  glmmTMB(ALA.LIN ~
            log(gonad.weight) +
            sex +
            TL + 
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchALA.LINmod1))

summary(perchALA.LINmod1)
# negative correlation with gonad weight
# positive correlation with total length

anova(perchALA.LINmod1,
      update(perchALA.LINmod1, ~-corral))  # corral significant

perchALA.LINpred <-
  ggemmeans(perchALA.LINmod1,
            terms = c("gonad.weight",
                      "TL")) %>% 
  rename(gonad.weight = x,
         TL = group)

png("Perch ALA-LIN Plot.png",
    width = 18,
    height= 10, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchALA.LINpred,
              aes(x = gonad.weight,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = TL),
              alpha = 0.3) +
  geom_line(data = perchALA.LINpred,
            aes(x = gonad.weight,
                y = predicted,
                colour = TL)) +
  geom_point(aes(x = gonad.weight,
                 y = ALA.LIN),
             size = 1) +
  scale_colour_viridis_d(option = "plasma",
                         name = "Total Length (cm)") +
  scale_fill_viridis_d(option = "plasma",
                       "Total Length (cm)") +
  labs(x = "Gonad Weight (g)",
       y = "ALA:LIN") +
  facet_wrap(~CCA1) +
  theme1

dev.off()

### PUFA/SFA+MUFA (feeding status) ----

plot(PUFA.SMUFA ~ corral, data = perch_FA2sex)

perchPUFA.SMUFAmod1 <- 
  glmmTMB(PUFA.SMUFA ~
            log(MPconcentration + 6) +
            TL +
            log(gonad.weight) +
            sex +
            (1 | corral),
          data =  perch_FA2sex)

plot(simulateResiduals(perchPUFA.SMUFAmod1))

summary(perchPUFA.SMUFAmod1)
# nothing is significant

anova(perchPUFA.SMUFAmod1,
      update(perchPUFA.SMUFAmod1, ~-corral))  # corral not significant

perchPUFA.SMUFApred <-
  ggemmeans(perchPUFA.SMUFAmod1,
            terms = c("gonad.weight")) %>% 
  rename(gonad.weight = x)

png("Perch PUFA-MUFA + SFA Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchPUFA.SMUFApred,
              aes(x = gonad.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchPUFA.SMUFApred,
            aes(x = gonad.weight,
                y = predicted)) +
  geom_point(aes(x = gonad.weight,
                 y = PUFA.SMUFA),
             size = 1) +
  labs(x = "Gonad Weight (g)",
       y = "PUFA/(SFA+MUFA)") +
  theme1

dev.off()

### ARA/LIN (n-6 comparison) ----

plot(ARA.LIN ~ corral, data = perch_FA2sex)

perchARA.LINmod1 <- 
  glmmTMB(ARA.LIN ~
            gonad.weight +
            sex +
            body.weight + 
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchARA.LINmod1))

summary(perchARA.LINmod1)
# nothing is significant

anova(perchARA.LINmod1,
      update(perchARA.LINmod1, ~-corral))  # corral not significant

perchARA.LINpred <-
  ggemmeans(perchARA.LINmod1,
            terms = c("gonad.weight")) %>% 
  rename(gonad.weight = x)

png("Perch ARA-LIN Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchARA.LINpred,
              aes(x = gonad.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchARA.LINpred,
            aes(x = gonad.weight,
                y = predicted)) +
  geom_point(aes(x = gonad.weight,
                 y = ARA.LIN),
             size = 1) +
  labs(x = "Gonad Weight (g)",
       y = "ARA:LIN") +
  theme1

dev.off()

### DHA/ALA (n-3 comparison) ----

plot(DHA.ALA ~ corral, data = perch_FA2sex)

perchDHA.ALAmod1 <- 
  glmmTMB(log(DHA.ALA) ~
            gonad.weight +
            sex +
            body.weight +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchDHA.ALAmod1))
# positive correlation with gonad weight

summary(perchDHA.ALAmod1)
# positive correlation with gonad weight
# weak negative correlation with body weight

anova(perchDHA.ALAmod1,
      update(perchDHA.ALAmod1, ~-corral))  # corral not significant

perchDHA.ALApred <-
  ggemmeans(perchDHA.ALAmod1,
            terms = c("gonad.weight")) %>% 
  rename(gonad.weight = x)

png("Perch DHA-ALA Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchDHA.ALApred,
              aes(x = gonad.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchDHA.ALApred,
            aes(x = gonad.weight,
                y = predicted)) +
  geom_point(aes(x = gonad.weight,
                 y = DHA.ALA),
             size = 1) +
  labs(x = "Gonad Weight (g)",
       y = "DHA:ALA") +
  theme1

dev.off()

### DHA/ARA (diet comparison) ----

plot(DHA.ARA ~ corral, data = perch_FA2sex)

perchDHA.ARAmod1 <- 
  glmmTMB(DHA.ARA ~
            TL +
            gonad.weight + 
            sex +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchDHA.ARAmod1))

summary(perchDHA.ARAmod1)
# nothing is significant

anova(perchDHA.ARAmod1,
      update(perchDHA.ARAmod1, ~-corral))  # corral not significant

perchDHA.ARApred <-
  ggemmeans(perchDHA.ARAmod1,
            terms = c("gonad.weight")) %>% 
  rename(gonad.weight = x)

png("Perch DHA-ARA Plot.png",
    width = 8.84,
    height= 5, 
    units = "cm",
    res = 500)

ggplot(perch_FA2sex) +
  geom_ribbon(data = perchDHA.ARApred,
              aes(x = gonad.weight,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = perchDHA.ARApred,
            aes(x = gonad.weight,
                y = predicted)) +
  geom_point(aes(x = gonad.weight,
                 y = DHA.ARA),
             size = 1) +
  labs(x = "Gonad Weight (g)",
       y = "DHA:ARA") +
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
                        log(MPconcentration + 6) + 
                        as.factor(date) +
                        (1 | corral),
                      data = zoop_FA_exp)

plotResiduals(simulateResiduals(zoopFAmod1))

summary(zoopFAmod1)

zoop_tot_FA_predict <- 
  ggemmeans(zoopFAmod1,
            terms = "date")

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

# Look at only mid and enpoints 

zoop_FA_exp <-
  zoop_FA %>% 
  filter(date != "2021-05-27")

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

zoopARAmod1 <- glmmTMB(log(C_20.4n.6) ~ 
                         log(MPconcentration + 6) +
                         as.factor(date) +
                          (1 | corral),
                        data = zoop_FA_exp)

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
                          log(MPconcentration + 6) +
                          as.factor(date) +
                          (1 | corral),
                        data = zoop_FA_exp)

plot(simulateResiduals(zoopEPAmod1))

summary(zoopEPAmod1)

zoopEPA_pred <- 
  ggemmeans(zoopEPAmod1,
            terms = c("date"),
            terms_to_colnames	= TRUE)

### DHA ----

zoopDHAmod1 <- glmmTMB(C_22.6n.3 ~ 
                          log(MPconcentration + 6) + 
                          as.factor(date) +
                          (1 | corral),
                        data = zoop_FA_exp)

plot(simulateResiduals(zoopDHAmod1))

summary(zoopDHAmod1)  # no effect

zoopDHA_pred <- ggemmeans(zoopDHAmod1,
                          terms = "date")

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

## PCA ----

# Pull out covariates

zoop_FA_conc_covariates <- zoop_FA[,c(1:3,51,52),]

zoop_FA_conc_totals <- zoop_FA[,c(17,26,35,45)]

summary(zoop_FA_conc_totals)

# Pull out fatty acid composition matrix

zoop_FA_conc_matrix <- zoop_FA[,c(7:16,18:25,27:34,36:44)]

### Run PCA ----
zoop_FA_conc_pca <- rda(zoop_FA_conc_matrix,
                         scale. = FALSE)

# Bar plot of relative eigenvalues
barplot(as.vector(zoop_FA_conc_pca$CA$eig)/sum(zoop_FA_conc_pca$CA$eig))

# Calculate percentage of variance explained by first 2 aaxes
sum((as.vector(zoop_FA_conc_pca$CA$eig)/sum(zoop_FA_conc_pca$CA$eig))[1:2])

summary(zoop_FA_conc_pca)

# 'Site' scores
zoop_FA_conc_PCA_site <- zoop_FA_conc_pca$CA$u

# 'Species' scores
zoop_FA_conc_PCA_species <- data.frame(zoop_FA_conc_pca$CA$v)
zoop_FA_conc_PCA_species$FA <- FA.names

zoop_FA_conc_PCA_site <- cbind(zoop_FA_conc_covariates,
                                zoop_FA_conc_PCA_site[, 1:2])

### Plot ----

png("Zooplankton FA Concentrations PCA.png",
    width = 19,
    height= 14, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_segment(data = zoop_FA_conc_PCA_species,
               aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_point(data = zoop_FA_conc_PCA_site,
             aes(x = PC1,
                 y = PC2,
                 colour = as.factor(MPconcentration),
                 shape = as.factor(date)),
             size = 4,
             alpha = 0.75) +
  geom_text(data = zoop_FA_conc_PCA_species,
            aes(x = PC1, 
                y = PC2, 
                label = FA),
            alpha = 0.5,
            size = 3,
            colour = "purple") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")")),
                      direction = -1) +
  scale_shape(name = "Date") +
  labs(x = "PC1",
       y = "PC2") +
  theme1

dev.off()

## nMDS ----

# Calculate distance matrix
zoop_FA_conc_diss <- as.matrix(vegdist(zoop_FA_conc_matrix, 
                                       method = "bray",
                                       na.rm = TRUE),
                               labels = TRUE)

NMDS_scree(zoop_FA_conc_diss)  # 4 dimensions looks good

set.seed(5465)

zoop_FA_conc_nMDS1 <- 
  metaMDS(zoop_FA_conc_diss,
          distance = "bray",
          k = 4,
          trymax = 250,
          wascores = TRUE,
          expand = TRUE,
          autotransform = FALSE)

# Shepards test/goodness of fit
goodness(zoop_FA_conc_nMDS1)
stressplot(zoop_FA_conc_nMDS1)

zoop_FA_conc_data.scores <- as.data.frame(scores(zoop_FA_conc_nMDS1))
zoop_FA_conc_data.scores2 <- cbind(zoop_FA_conc_data.scores,
                              zoop_FA_conc_covariates)

zoop_FA_conc_scores <- `sppscores<-`(zoop_FA_conc_nMDS1, zoop_FA_conc_matrix)

zoop_FA_conc_variable_scores <- 
  as.data.frame(zoop_FA_conc_scores$species)

zoop_FA_conc_variable_scores$FA <- FA.names

### Generate hulls ----

zoop_FA_conc_data.scores2$MPconcentration <- 
  as.factor(zoop_FA_conc_data.scores2$MPconcentration)

zoop_FA_conc_data.scores2$date <- 
  as.factor(zoop_FA_conc_data.scores2$date)

zoop_FA_conc_hulls <- data.frame()

for(i in 1:length(unique(zoop_FA_conc_data.scores2$MPconcentration))) {
  hull <-
    zoop_FA_conc_data.scores2[zoop_FA_conc_data.scores2$MPconcentration ==
                           unique(zoop_FA_conc_data.scores2$MPconcentration)[i],
    ][chull(zoop_FA_conc_data.scores2[zoop_FA_conc_data.scores2$MPconcentration ==
                                   unique(zoop_FA_conc_data.scores2$MPconcentration)[i],
                                 c(1:2)]),]
  zoop_FA_conc_hulls <- rbind(zoop_FA_conc_hulls, hull)
}

### Plot ----

png("Zooplankton FA Concentrations MDS Plot.png",
    width = 19,
    height= 12, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = zoop_FA_conc_hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   fill = MPconcentration,
                   colour = MPconcentration),
               alpha = 0.3,
               size = 0.5) +
  geom_segment(data = zoop_FA_conc_variable_scores,
               aes(x = 0, y = 0, xend = MDS1*0.9, yend = MDS2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = zoop_FA_conc_hulls,
             aes(x = NMDS1,
                 y = NMDS2,
                 shape = date,
                 colour = MPconcentration),
             alpha = 0.5,
             size = 3) +
  geom_text(data = zoop_FA_conc_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = FA),
            alpha = 0.9,
            size = 8 / .pt,
            colour = "purple") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      direction = -1,
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_fill_brewer(type = "div",
                    palette = "RdYlGn",
                    direction = -1,
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_shape(name = "Date") +
  theme1

dev.off()



# Fish food ----

png("Fish Food Essential FAs Concentrations Plot.png",
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

