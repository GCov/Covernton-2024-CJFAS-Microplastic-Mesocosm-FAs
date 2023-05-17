library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)

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

#### Compare body size with MP concentration ----

fish2021.mod1 <-
  aov(body.weight ~ corral,
      data = perch2021)

summary(fish2021.mod1)

plot(residuals(fish2021.mod1) ~ fitted(fish2021.mod1))
abline(0,0)

# Corral is sig (p<0.001)

TukeyHSD(fish2021.mod1)

## Differences in H-B, H-C, H-D, 

fish2021.pred <- 
  data.frame(corral = levels(perch2021$corral))

fish2021.pred$fit <- predict(fish2021.mod1,
                             newdata = fish2021.pred,
                             re.form = NA)

fish2021.pred <- left_join(fish2021.pred,
                           pop,
                           by = "corral")

fish2021.labs <- 
  perch2021 %>% 
  group_by(corral, MPconcentration) %>% 
  summarize(y = max(na.omit(body.weight)))

fish2021.labs$lab <-
  c("a", "a", "a", "ab", "ab", "ab", "b", "ab")

png("2021 Perch Weights.png",
    width = 8.84,
    height= 6, 
    units = "cm",
    res = 600)

ggplot(perch2021) +
  geom_boxplot(aes(x = MPconcentration,
                   y = body.weight,
                   fill = corral),
               alpha = 0.75) +
  geom_text(data = fish2021.labs,
            aes(x = MPconcentration,
                y = y + 0.5,
                label = lab,
                group = corral),
            size = 10/.pt,
            position = position_dodge(width = 1)) +
  labs(x = expression(paste("MP Exposure Concentration (particles"~L^-1*")")),
       y = "Body Weight (g)") +
  scale_colour_viridis_d(option = "turbo",
                         name = "Corral") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Corral") +
  scale_x_continuous(trans = "log1p",
                     breaks = sort(unique(perch2021$MPconcentration))) +
  theme1

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
