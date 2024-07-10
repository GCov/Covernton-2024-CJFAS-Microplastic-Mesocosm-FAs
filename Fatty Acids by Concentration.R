# Setup ----

# Load libraries

library(ggplot2)  # for plotting
library(tidyr)  # for data wrangling
library(dplyr)  # for data wrangling
library(glmmTMB)  # for LMs
library(DHARMa)  # for assessing LMs
library(emmeans)  # for interpreting LMs
library(ggeffects)  # intaerpreting LMs
library(vegan)  # for multivariate modeling
library(ggrepel)  # for text plotting
library(cowplot)  # for making plot panels

# Set a universal ggplot theme to use for all plots

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

FA.names <-
  c(
    "12:0",
    "13:0",
    "14:0",
    "15:0",
    "16:0",
    "17:0",
    "18:0",
    "20:0",
    "22:0",
    "24:0",
    "12:1",
    "14:1",
    "16:1n-7",
    "16:1n-9",
    "18:1n-7",
    "18:1n-9",
    "20:1n-9",
    "22:1n-9",
    "18:2n-6",
    "18:3n-6",
    "20:2n-6",
    "20:3n-6",
    "20:4n-6",
    "22:2n-6",
    "22:4n-6",
    "22:5n-6",
    "18:3n-3",
    "18:4n-3",
    "20:3n-3",
    "20:4n-3",
    "20:5n-3",
    "22:5n-3",
    "22:6n-3",
    "24:5n-3",
    "24:6n-3"
  )

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

# Add population data ----

fish_pop <- read.csv("fish_pop.csv", header = TRUE)

str(fish_pop)

fish_pop$corral <- as.factor(fish_pop$corral)
fish_pop$MPconcentration <- as.numeric(fish_pop$MPconcentration)

fish_pop$perch.surv <- with(fish_pop, YP.end / YP.start)

# Combine with FA data

perch_FA2 <- left_join(perch_FA2, fish_pop,
                       by = c("corral", "MPconcentration"))

# Perch Analysis----

# Reduce dataset to fish with sex and gonad weight
# For later analyses

perch_FA2sex <-
  perch_FA2 %>%
  filter(!is.na(sex) &
           !is.na(gonad.weight))

nrow(perch_FA2sex)  # 18 fish remaining

## Total FAs ----

perchFAmod1 <- glmmTMB(total_FAs ~
                         body.weight +
                         log(MPconcentration + 6) +
                         (1 | corral),
                       data = perch_FA2)

plot(simulateResiduals(perchFAmod1))

summary(perchFAmod1)
# body weight marginally significant

anova(perchFAmod1, update(perchFAmod1, ~ -corral))  # corral not significant

## Relationship of all FAs with total FAs ----

### Put data into long form

perch_FAlong <-
  perch_FA2 %>%
  pivot_longer(cols = c(7:16, 18:25, 27:34, 36:44),
               names_to = "FAs")

perch_FAlong$FAs <- as.factor(perch_FAlong$FAs)

### Exploratory plots ----

# total FAs
ggplot(perch_FAlong) +
  geom_smooth(aes(x = total_FAs,
                  y = value,
                  colour = FAs),
              method = "lm") +
  geom_point(aes(x = total_FAs,
                 y = value,
                 colour = FAs)) +
  facet_wrap( ~ FAs) +
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
  facet_wrap( ~ FAs) +
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
  facet_wrap( ~ FAs) +
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
  facet_wrap( ~ FAs) +
  theme1

## Individual FAs -----

### Treatment effect on HUFAs ----

#### EPA ----

perchEPAmod1 <-
  glmmTMB(C_20.5n.3 ~
            log(MPconcentration + 6) +
            body.weight +
            (1 | corral),
          data = perch_FA2)

plot(simulateResiduals(perchEPAmod1))

summary(perchEPAmod1)
# positively correlated with body weight

perchEPA_pred <-
  ggemmeans(perchEPAmod1,
            terms = c("body.weight")) %>%
  rename(body.weight = x)

png(
  "Perch EPA Plot.png",
  width = 8.4,
  height = 5,
  units = "cm",
  res = 500
)

ggplot(perch_FA2) +
  geom_ribbon(
    data = perchEPA_pred,
    aes(x = body.weight,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.3
  ) +
  geom_line(data = perchEPA_pred,
            aes(x = body.weight,
                y = predicted)) +
  geom_point(aes(x = body.weight,
                 y = C_20.5n.3),
             size = 1) +
  scale_colour_viridis_d(option = "plasma",
                         name = "Total Length (cm)") +
  labs(x = "Body Weight (g)",
       y = expression(paste("EPA Concentration (mg " ~ g ^ -1 * ")"))) +
  theme1

dev.off()

#### DHA ----

perchDHAmod1 <-
  glmmTMB(C_22.6n.3 ~
            log(MPconcentration + 6) +
            body.weight +
            (1 | corral),
          data = perch_FA2)

plot(simulateResiduals(perchDHAmod1))

summary(perchDHAmod1)
# positively correlated with MP concentration

perchDHA_pred <-
  ggemmeans(perchDHAmod1,
            terms = c("MPconcentration [0.001:29240, by = 200]")) %>%
  rename(MPconcentration = x)

png(
  "Perch DHA Plot.png",
  width = 8.84,
  height = 5,
  units = "cm",
  res = 300
)

ggplot(perch_FA2) +
  geom_ribbon(
    data = perchDHA_pred,
    aes(x = MPconcentration,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.3
  ) +
  geom_line(data = perchDHA_pred,
            aes(x = MPconcentration,
                y = predicted)) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3),
             size = 1) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  labs(x = expression(paste("Exposure Concentration (MPs" ~
                              L ^ -1 * ")")),
  y = expression(paste("Exposure Concentration (MPs" ~
                         L ^ -1 * ")"))) +
  theme1

dev.off()

#### ARA ----

perchARAmod1 <-
  glmmTMB(C_20.4n.6 ~
            body.weight +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA2)

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
  geom_ribbon(
    data = perchALA_pred,
    aes(
      x = gonad.weight,
      ymin = conf.low,
      ymax = conf.high,
      fill = body.weight
    ),
    alpha = 0.3
  ) +
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
       y = expression(paste("ALA Concentration (mg " ~ g ^ -1 * ")"))) +
  theme1

#### LA ----

perchLAmod1 <-
  glmmTMB(C_18.2n.6 ~
            log(gonad.weight) +
            body.weight +
            sex +
            (1 | corral),
          data = perch_FA2sex)

plot(simulateResiduals(perchLAmod1))

summary(perchLAmod1)
# weak negative correlation with gonad weight

perchLA_pred <-
  ggemmeans(perchLAmod1,
            terms = c("gonad.weight")) %>%
  rename(gonad.weight = x)

LAplot <-
  ggplot(perch_FA2sex) +
  geom_ribbon(
    data = perchLA_pred,
    aes(x = gonad.weight,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.3
  ) +
  geom_line(data = perchLA_pred,
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
       y = expression(paste("LA Concentration (mg " ~ g ^ -1 * ")"))) +
  theme1

#### ALA and LA plot ----

png(
  "Perch ALA and LA Plot.png",
  width = 18,
  height = 6,
  units = "cm",
  res = 300
)

plot_grid(
  ALAplot,
  LAplot,
  ncol = 2,
  align = "v",
  axis = "r",
  labels = c("A", "B")
)

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

png(
  "Perch DHA Plot 2.png",
  width = 8.84,
  height = 5,
  units = "cm",
  res = 300
)

ggplot(perch_FA2sex) +
  geom_ribbon(
    data = perchDHA_pred2,
    aes(x = body.weight,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.3
  ) +
  geom_line(data = perchDHA_pred2,
            aes(x = body.weight,
                y = predicted)) +
  geom_point(aes(x = body.weight,
                 y = C_22.6n.3),
             size = 1) +
  labs(x = "Body Weight (g)",
       y = expression(paste("DHA Concentration (mg " ~ g ^ -1 * ")"))) +
  theme1

dev.off()

# Zooplankton Analysis----

## Plot different FAs

names(zoop_FA)

plot(total_SFAs ~ date, data = zoop_FA)
plot(total_MUFAs ~ date, data = zoop_FA)
plot(total_N.3_PUFA ~ date, data = zoop_FA)
plot(total_N.6_PUFA ~ date, data = zoop_FA)

# Separate mid and endpoints

zoop_FA_exp <-
  zoop_FA %>%
  filter(date != "2021-05-27")

## Total FAs----

zoopFAmod1 <- glmmTMB(total_FAs ~
                        log(MPconcentration + 6) +
                        as.factor(date) +
                        (1 | corral),
                      data = zoop_FA_exp)

plotResiduals(simulateResiduals(zoopFAmod1))

summary(zoopFAmod1)


## Individual FAs ----

### ARA ----

zoopARAmod1 <- glmmTMB(log(C_20.4n.6) ~
                         log(MPconcentration + 6) +
                         as.factor(date) +
                         (1 | corral),
                       data = zoop_FA_exp)

plot(simulateResiduals(zoopARAmod1))

summary(zoopARAmod1)  # no effect

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

zoopEPA_pred

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

zoopDHA_pred


## Plot compositon by absolute concentration ----

zoop_FA_long <-
  zoop_FA %>%
  pivot_longer(cols = c(17, 26, 35, 45),
               names_to = "FA")

zoop_FA_long$FA <- as.factor(zoop_FA_long$FA)

levels(zoop_FA_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")


zoop_FA_long$timepoint <-
  as.factor(zoop_FA_long$date)

levels(zoop_FA_long$timepoint) <-
  c("Day -6",
    "Day 34",
    "Day 68")

png(
  "Zooplankton Absolute Composition Plot.png",
  width = 18,
  height = 12,
  units = "cm",
  res = 600
)

ggplot(zoop_FA_long) +
  geom_col(aes(
    x = reorder(ID, as.numeric(date)),
    y = value,
    fill = FA,
    colour = timepoint
  ),
  size = 0.75) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  scale_colour_viridis_d(name = "Time Point",
                         option = "plasma") +
  labs(x = "Sample",
       y = expression(paste("Concentration (mg " ~ g ^ -1 * ")"))) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap( ~ MPconcentration,
              scales = "free_x") +
  theme1 +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

dev.off()

## Plot by Individual FA ----

# Put data into long form

str(zoop_FA)

zoop_FA_long <-
  zoop_FA %>% 
  pivot_longer(cols = c(7:16, 18:25, 27:34, 36:44),
               names_to = "FA") %>% 
  group_by(FA) %>% 
  filter(mean(value)/mean(total_FAs) > 0.01) %>% 
  ungroup() %>% 
  mutate(time = as.factor(as.character(date)))

levels(zoop_FA_long$time) <- c("Day -6",
                               "Day 34",
                               "Day 68")

zoop_FA_names <- 
  c("14:0",
    "16:0",
    "16:1n-7",
    "18:0",
    "18:1n-7",
    "18:1n-9",
    "18:2n-6",
    "18:3n-3",
    "18:3n-6",
    "18:4n-3",
    "20:4n-6",
    "20:5n-3",
    "22:1n-9",
    "22:5n-6",
    "22:6n-3")
  
zoop_FA_long$FA <- as.factor(zoop_FA_long$FA)

levels(zoop_FA_long$FA) <- zoop_FA_names

png(
  "Zoop Individual FAs.png",
  width = 18,
  height = 18,
  units = "cm",
  res = 300
)

ggplot(zoop_FA_long) +
  geom_col(aes(x = corral,
               y = value,
               fill = FA)) +
  facet_grid(time~MPconcentration, 
             scales = "free_x", 
             space = "free") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis_d(option = "turbo",
                       name = "Fatty Acid") +
  labs(y = expression(paste("Concentration (mg " ~ g ^ -1 * ")")),
       x = "") +
  theme1 +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

dev.off()



# Fish food ----

png(
  "Fish Food Essential FAs Concentrations Plot.png",
  width = 9,
  height = 9,
  units = "cm",
  res = 600
)

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
       y = expression(paste("Concentration (mg " ~ g ^ -1 * ")"))) +
  scale_fill_viridis_d(name = "") +
  theme1 +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

dev.off()

