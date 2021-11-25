library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(MuMIn)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"), 
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

## Load data

perch_diet <- read.csv("perch_diet.csv", header = TRUE)

perch_diet$dose <- as.numeric(perch_diet$dose)

## Plot totals

ggplot(data = perch_diet,
       aes(x = dose,
           y = total.animals)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

## Put data into long format

perch_diet_long <- 
  perch_diet %>%
  pivot_longer(names(perch_diet[12:17]),
               names_to = "taxa",
               values_to = "count")

str(perch_diet_long)

perch_diet_long$taxa <- as.factor(perch_diet_long$taxa)
perch_diet_long$ID <- as.factor(perch_diet_long$ID)
perch_diet_long$count <- as.numeric(perch_diet_long$count)

levels(perch_diet_long$taxa) <-
  c("Adult Insects",
    "Amphipods",
    "Larval Chironomids",
    "Cladocerans",
    "Cyclopoid Copepods",
    "Larval Odonates")

perch_diet_long$treatment <-
  as.factor(perch_diet_long$corral)

levels(perch_diet_long$treatment) <-
  c("0(B)", "414", "29,240", "100", "6", "7,071", "0(H)", "1,710")

## Plot by taxa

png("Perch Diet Plot by Taxa.png",
    width = 25,
    height= 12, 
    units = "cm",
    res = 600)

ggplot(perch_diet_long) +
  geom_col(aes(x = ID,
               y = count,
               fill = reorder(taxa, 1/(count+1), mean)),
           colour = "black",
           size = 0.25) +
  labs(x = expression(paste("Dose (MPs"~L^-1*")")),
       y = "Number of Individuals (log scale)") +
  scale_fill_brewer(type = "qual",
                    name = "Taxa",
                    palette = "Set1") +
  scale_y_continuous(expand = c(0.005,0.005),
                     breaks = c(0, 1, 10, 1000)) +
  coord_trans(y = "log1p") +
  facet_grid(.~reorder(treatment, dose, mean),
             scales = "free_x",
             switch = "x",
             space = "free_x") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

dev.off()

## Try a glmm

# Poisson
mod1 <- glmmTMB(total.animals ~ log(dose + 1) + body.length +  (1 | corral),
                family = poisson(link = "log"),
                data = perch_diet)  
summary(mod1)
res1 <- simulateResiduals(mod1)
plot(res1)

# NB with linear variance
mod2 <- glmmTMB(total.animals ~ log(dose + 1) + body.length +  (1 | corral),
                family = nbinom1(link = "log"),
                data = perch_diet)
summary(mod2)
res2 <- simulateResiduals(mod2)
plot(res2)

# NB with quadratic variance
mod3 <- glmmTMB(total.animals ~ log(dose + 1) + body.length +  (1 | corral),
                family = nbinom2(link = "log"),
                data = perch_diet)
summary(mod3)
res3 <- simulateResiduals(mod3)
plot(res3)

AICc(mod2, mod3)  # mod3 is the better fit

## Plot model predictions

simdata <- data.frame(dose = seq(0, max(perch_diet$dose), length.out = 1000),
                      body.length = rep(mean(perch_diet$body.length), 1000),
                      corral = rep(NA, 1000))

simdata$mean <- predict(mod3, simdata)

simdata$se <- as.numeric(predict(mod3,simdata, se.fit = TRUE)$se.fit)

simdata$upper <- with(simdata, mean + 1.96*se)
simdata$lower <- with(simdata, mean - 1.96*se)

ggplot() +
  geom_ribbon(data = simdata,
              aes(x = dose,
                  ymin = lower,
                  ymax = upper),
              fill = "blue",
              alpha = 0.5) +
  geom_line(data = simdata,
            aes(x = dose,
                y = mean),
            size = 0.5) +
  geom_point(data = perch_diet,
             aes(x = dose,
                 y = total.animals)) +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1
