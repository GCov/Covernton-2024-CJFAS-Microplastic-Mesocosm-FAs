library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 9,
                        family = "serif"),
    axis.text = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing = unit(0, "cm")
  )

#### Load data ----

biometrics <- read.csv("perch_biometrics.csv",
                       header = TRUE,
                       stringsAsFactors = TRUE)

pop <- read.csv("fish_pop.csv",
                header = TRUE,
                stringsAsFactors = TRUE)

perch2021 <- left_join(biometrics, pop, by = "corral")

# Remove those without data (mostly non-perch fish)

perch2021 <- 
  perch2021 %>%
  filter(!is.na(TL))

# Summarize

perch2021 %>% 
  group_by(MPconcentration, YP.end) %>% 
  summarize(mean.body.weight = mean(body.weight),
            sd.body.weight = sd(body.weight),
            mean.TL = mean(TL),
            sd.TL = sd(TL))

#### Compare body size with MP concentration ----

fish2021.mod1 <-
  aov(body.weight ~ corral,
      data = perch2021)

summary(fish2021.mod1)

plot(residuals(fish2021.mod1) ~ fitted(fish2021.mod1))
abline(0,0)

# Corral is sig (p<0.001)

TukeyHSD(fish2021.mod1)

## Differences in H-B, H-D, H-I

fish2021.pred <- 
  data.frame(corral = levels(perch2021$corral))

fish2021.pred$fit <- predict(fish2021.mod1,
                             newdata = fish2021.pred,
                             re.form = NA)

fish2021.pred <- left_join(fish2021.pred,
                           pop,
                           by = "corral")

perch2021$treatment <-
  perch2021$MPconcentration

perch2021$treatment[perch2021$treatment == "0"] <- 
  with(perch2021[perch2021$treatment == "0",],
       ifelse(corral == "B",
              "0 (1)",
              "0 (2)")) 

perch2021$treatment <- as.factor(perch2021$treatment)

fish2021.labs <- 
  perch2021 %>% 
  group_by(corral, MPconcentration, treatment) %>% 
  summarize(y = max(na.omit(body.weight)))

fish2021.labs$lab <-
  c("a", "ab", "a", "ab", "ab", "ab", "b", "a")

png("2021 Perch Weights.png",
    width = 8.84,
    height= 6, 
    units = "cm",
    res = 600)

ggplot(perch2021) +
  geom_boxplot(aes(x = MPconcentration,
                   y = body.weight,
                   fill = reorder(treatment,
                                  MPconcentration,
                                  mean)),
               alpha = 0.75) +
  geom_text(data = fish2021.labs,
            aes(x = MPconcentration,
                y = y + 0.5,
                label = lab,
                group = reorder(treatment,
                                MPconcentration,
                                mean)),
            size = 10/.pt,
            position = position_dodge(width = 1)) +
  labs(x = expression(paste("MP Exposure Concentration (particles"~L^-1*")")),
       y = "Body Weight (g)") +
  scale_fill_viridis_d(option = "plasma",
                       name = "Treatment") +
  scale_x_continuous(trans = "log1p",
                     breaks = sort(unique(perch2021$MPconcentration))) +
  theme1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

dev.off()

ggplot(perch2021) +
  geom_density(aes(x = FL,
                   fill = corral),
               alpha = 0.5) +
  facet_grid(MPconcentration ~ .) +
  labs(x = "Fork Length (cm)",
       y = "Density") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Corral") +
  scale_x_continuous(limits = c(6, 11.5),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme1

perch2021$FultonK <- with(perch2021,
                          body.weight / (TL ^ 3))

fish2022.mod3 <-
  glmmTMB(FultonK ~ corral,
          data = perch2021)

summary(fish2022.mod3)

plot(simulateResiduals(fish2022.mod3))

ggplot(perch2021) +
  geom_density(aes(x = FultonK,
                   fill = corral),
               alpha = 0.5) +
  facet_grid(MPconcentration ~ .) +
  labs(x = "Condition Factor (K)",
       y = "Density") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Corral") +
  scale_x_continuous(limits = c(0.0075, 0.0125),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme1


## Neither fork length nor condition affected by MP concentration

## Gonad weight

fish2021.mod2 <-
  aov(gonad.weight ~ corral,
      data = perch2021)

summary(fish2021.mod2)

plot(residuals(fish2021.mod1) ~ fitted(fish2021.mod1))
abline(0,0)

## Check whether regression is better

names(perch2021)

fish2021.mod2 <-
  glmmTMB(body.weight ~ log(MPconcentration + 6) + YP.end,
      data = perch2021)

plot(simulateResiduals(fish2021.mod2))

summary(fish2021.mod2)

plot(body.weight ~ log(MPconcentration + 6), data = perch2021)
lines(predict(fish2021.mod2, re.form = NA) ~ log(perch2021$MPconcentration + 6))

plot(residuals(fish2021.mod2) ~ predict(fish2021.mod2))

plotResiduals(fish2021.mod2, perch2021$date.collected)
plotResiduals(fish2021.mod2, perch2021$YP.end)

perch.body.predict <-
  as.data.frame(ggemmeans(fish2021.mod2,
                          c("MPconcentration [0.001:29240,by = 1000]", 
                            "YP.end"))) %>% 
  rename(MPconcentration = x,
         YP.end = group)

ggemmeans(fish2021.mod2,
          c("YP.end"))

perch.body.predict$YP.end <- as.character(perch.body.predict$YP.end)
perch.body.predict$YP.end <- as.numeric(perch.body.predict$YP.end)

png("2021 Perch Weights 2.png",
    width = 8.84,
    height= 8, 
    units = "cm",
    res = 600)

ggplot(perch2021) +
  geom_ribbon(data = perch.body.predict,
              aes(x = MPconcentration,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = YP.end,
                  group = as.factor(YP.end)),
              alpha = 0.4) +
  geom_line(data = perch.body.predict,
            aes(x = MPconcentration,
                y = predicted,
                colour = YP.end,
                group = as.factor(YP.end))) +
  geom_point(aes(x = MPconcentration,
                 y = body.weight,
                 fill = YP.end),
             shape = 21) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 6, 100, 414, 1710, 7071, 29240)) +
  scale_fill_viridis_c(option = "plasma",
                       name = "Surviving Yellow Perch") +
  scale_colour_viridis_c(option = "plasma",
                         name = "Surviving Yellow Perch") +
  labs(x = expression(paste("Microplastic exposure concentration (particles"~L^-1*")")),
       y = "Body Weight (g)") +
  theme1 +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2))

dev.off()
