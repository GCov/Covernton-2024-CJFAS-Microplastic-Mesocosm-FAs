# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)

# Load data ----

FAs <- read.csv("FAs.csv", header = TRUE)
str(FAs)
FAs$ID <- as.factor(FAs$ID)
FAs$corral <- as.factor(FAs$corral)
FAs$sample.type <- as.factor(FAs$sample.type)

FAs <- left_join(FAs, perch_diet[,c(2:3)])

# Separate out by sample type ----

perch_FA <- subset(FAs, sample.type == "perch")
zoop_FA <- subset(FAs, sample.type == "zooplankton")
food_FA <- subset(FAs, sample.type == "fish food")

# Make plots by total FAs and DHA and EPA ----

## Perch plot ----

ggplot(perch_FA) +
  geom_point(aes(x = MPconcentration,
                 y = Total_FAs)) +
  geom_smooth(aes(x = MPconcentration,
                  y = Total_FAs),
              colour = "black",
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = EPA.DHA),
             colour = "red") +
  geom_smooth(aes(x = MPconcentration,
                  y = EPA.DHA),
              colour = "red",
              method = "lm") +
  labs(x = expression(paste("MP exposure cocentration (particles"~L^-1*")")),
       y = expression(paste("Fatty acids concentration (mg"*g^-1*")"))) +
  theme1

## Zoop plot ----

ggplot(zoop_FA) +
  geom_point(aes(x = MPconcentration,
                 y = Total_FAs)) +
  geom_smooth(aes(x = MPconcentration,
                  y = Total_FAs),
              colour = "black",
              method = "lm") +
  geom_point(aes(x = MPconcentration,
                 y = EPA.DHA),
             colour = "red") +
  geom_smooth(aes(x = MPconcentration,
                  y = EPA.DHA),
              colour = "red",
              method = "lm") +
  labs(x = expression(paste("MP exposure cocentration (particles"~L^-1*")")),
       y = expression(paste("Fatty acids concentration (mg"*g^-1*")"))) +
  theme1
ggplot(zoop_FA) +
  geom_boxplot(aes(x = corral,
                   y = Total_FAs)) +
  geom_boxplot(aes(x = corral,
                   y = EPA.DHA),
               fill = "red") +
  labs(x = "Corral",
       y = expression(paste("Fatty acids concentration (mg"*g^-1*")"))) +
  theme1
