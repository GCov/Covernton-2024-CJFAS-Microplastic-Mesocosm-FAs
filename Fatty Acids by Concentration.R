# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)

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

# Separate out by sample type ----

perch_FA <- subset(FAs, sample.type == "perch")
zoop_FA <- subset(FAs, sample.type == "zooplankton")
food_FA <- subset(FAs, sample.type == "fish food")

# Assign corral a MP concentration of 0 for zoop_FA

zoop_FA[zoop_FA$corral == "A",]$MPconcentration <- 24

# Add in perch biometrics data ----

perch_biometrics <- read.csv("perch_biometrics.csv", header = TRUE)

str(perch_biometrics)

perch_biometrics$corral <- as.factor(perch_biometrics$corral)
perch_biometrics$ID <- as.factor(perch_biometrics$ID)
perch_biometrics$sex <- as.factor(perch_biometrics$sex)

# Combine with perch FA data

perch_FA2 <- left_join(perch_FA, perch_biometrics, 
                       by = c("ID", "corral"))

# Analysis and plots for total FAs----

## Perch ----

# Convert date to a centered numeric value

perch_FA2$date2 <- scale(as.numeric(perch_FA2$date), center = TRUE)

perchFAmod1 <- glmmTMB(total_FAs ~ 
                         scale(MPconcentration, center = TRUE) +
                         (1 | corral),
                       data = perch_FA2)

plotResiduals(simulateResiduals(perchFAmod1))

summary(perchFAmod1)

perch_tot_FA_predict <- predict(perchFAmod1,
                                re.form = NA, 
                                se.fit = TRUE)

perch_tot_FA_predict$upper <- with(perch_tot_FA_predict,
                                   fit + 1.98 * se.fit)

perch_tot_FA_predict$lower <- with(perch_tot_FA_predict,
                                   fit - 1.98 * se.fit)

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
  labs(x = expression(paste("MP exposure cocentration (particles"~L^-1*")")),
       y = expression(paste("Total fatty acids (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_brewer(palette = "Set1",
                      name = "Corral") +
  theme1

## Zooplankton ----

# Convert date to a factor

zoop_FA$date2 <- as.factor(zoop_FA$date)

zoopFAmod1 <- glmmTMB(total_FAs ~
                        MPconcentration*date2 +
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
            size = 2) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Total fatty acids (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_brewer(palette = "Set1",
                      name = "Corral") +
  facet_wrap(~date) +
  theme1

# Analysis and plots according to DHA vs. EPA -----

## Perch ----

perchEPAmod1 <- glmmTMB(C_20.5n.3 ~ 
                          log(MPconcentration + 1) +
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchEPAmod1))

summary(perchEPAmod1)

perchDHAmod1 <- glmmTMB(C_22.6n.3 ~ 
                          log(MPconcentration + 1) + 
                          date2 +
                          (1 | corral),
                        data = perch_FA2)

plotResiduals(simulateResiduals(perchDHAmod1))

summary(perchDHAmod1)

colours <- c("DHA" = "green",
             "EPA" = "blue")

ggplot(perch_FA) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.5n.3,
                 colour = "EPA")) +
  geom_line(aes(x = MPconcentration,
                y = predict(perchEPAmod1),
                colour = "EPA")) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3,
                 colour = "DHA")) +
  geom_line(aes(x = MPconcentration,
                y = predict(perchDHAmod1),
                colour = "DHA")) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = colours,
                      name = "") +
  theme1

## Zooplankton ----

ggplot(zoop_FA) +
  geom_point(aes(x = MPconcentration,
                 y = C_20.5n.3,
                 colour = "EPA")) +
  geom_smooth(aes(x = MPconcentration,
                  y = C_20.5n.3,
                  colour = "EPA"),
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3,
                 colour = "DHA")) +
  geom_smooth(aes(x = MPconcentration,
                  y = C_22.6n.3,
                  colour = "DHA"),
              method = "lm") +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = colours,
                      name = "") +
  facet_grid(~date) +
  theme1

## Fish food ----

ggplot(food_FA) +
  geom_boxplot(aes(x = 1,
                   y = C_20.5n.3,
                   fill = "EPA")) +
  geom_boxplot(aes(x = 2,
                   y = C_22.6n.3,
                   fill = "DHA")) +
  labs(x = "",
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_fill_manual(values = colours,
                      name = "") +
  theme1 +
  theme(axis.text.x = element_blank())

