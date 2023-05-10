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
    text = element_text(size = 7,
                        family = "serif"),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.key.size = unit(0.5, "cm")
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
  glmmTMB(body.weight ~ corral,
          data = perch2021)

summary(fish2021.mod1)

plot(simulateResiduals(fish2021.mod1))

fish2021.mod2 <-
  glmmTMB(body.weight ~ 1,
          data = perch2021)

anova(fish2021.mod1,
      fish2021.mod2)  # Corral is sig (p<0.001)

plot(emmeans(fish2021.mod1,
             specs = "corral"))

## Corral H (one of the controls) had smaller fish than all other corrals,
## except possibly E(100).
## Corral B (the other control) had larger fish than E(100)

fish2021.pred <- 
  data.frame(corral = levels(perch2021$corral))

fish2021.pred$fit <- predict(fish2021.mod1,
                             newdata = fish2021.pred,
                             re.form = NA)

fish2021.pred <- left_join(fish2021.pred,
                           pop,
                           by = "corral")

png("2021 Perch Weights.png",
    width = 12,
    height= 12, 
    units = "cm",
    res = 500)

ggplot(perch2021) +
  geom_density(aes(x = body.weight,
                   fill = corral),
               alpha = 0.5) +
  geom_vline(data = fish2021.pred,
             aes(xintercept = fit,
                 colour = corral),
             linetype = "dashed") +
  facet_grid(MPconcentration ~ .) +
  labs(x = "Body Weight (g)",
       y = "Density") +
  scale_colour_viridis_d(option = "turbo",
                         name = "Corral") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Corral") +
  scale_x_continuous(limits = c(1, 16), 
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
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