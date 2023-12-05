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

# Set a universal ggplot theme to use for all plots

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


## Load data ----

perch_diet <- read.csv("perch_diet.csv", header = TRUE)

perch_diet$MPconcentration <- as.numeric(perch_diet$MPconcentration)
perch_diet$corral <- as.factor(perch_diet$corral)

# Add number of surviving perch in each corral
# *Note that this relies on running the Perch Survival and Growth script first

perch_diet <-
  left_join(perch_diet,
            perch2021[, c(2, 13)])

# Individual count plots ----

## Plot totals

# Everything

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = total.animals)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles" ~ L ^ -1 * ")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

# Cladocerans

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = cladocera)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles" ~ L ^ -1 * ")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

# Copepods

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = cyclopoida)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles" ~ L ^ -1 * ")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

# Pupal chironomids

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = chironomid.pupa)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles" ~ L ^ -1 * ")")),
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
  c(
    "Amphipods",
    "Chironomid Larvae",
    "Chironomid Pupae",
    "Cladocerans",
    "Cyclopoid Copepods",
    "Larval Odonates"
  )

perch_diet_long$treatment <-
  as.factor(perch_diet_long$corral)

levels(perch_diet_long$treatment) <-
  c("0(1)", "414", "29,240", "100", "6", "7,071", "0(2)", "1,710")

## Plot by taxa ----

png(
  "Perch Diet Plot by Taxa.png",
  width = 18,
  height = 8,
  units = "cm",
  res = 600
)

ggplot(perch_diet_long) +
  geom_col(aes(
    x = ID,
    y = count,
    fill = reorder(taxa, 1 / (count + 1), mean)
  ),
  colour = "black",
  size = 0.25) +
  labs(x = expression(paste(
    "MP exposure concentration (particles" ~ L ^ -1 * ")"
  )),
  y = "Number of Individuals (log scale)") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Taxa") +
  scale_y_continuous(expand = c(0.005, 0.005),
                     breaks = c(0, 1, 10, 100, 1000)) +
  coord_trans(y = "log1p") +
  facet_grid(
    . ~ reorder(treatment, MPconcentration, mean),
    scales = "free_x",
    switch = "x",
    space = "free_x"
  ) +
  theme1 +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm")
  )

dev.off()


## Plot by relative abundance ----

perch_relabund <- perch_diet
perch_relabund[, c(12:17)] <-
  decostand(perch_relabund[, c(12:17)],
            method = "total")

perch_relabund_long <-
  perch_relabund %>%
  pivot_longer(names(perch_relabund[12:17]),
               names_to = "taxa",
               values_to = "count")

perch_relabund_long$taxa <- as.factor(perch_relabund_long$taxa)
perch_relabund_long$ID <- as.factor(perch_relabund_long$ID)
perch_relabund_long$count <- as.numeric(perch_relabund_long$count)

levels(perch_relabund_long$taxa) <-
  c(
    "Amphipods",
    "Chironomid Larvae",
    "Chironomid Pupae",
    "Cladocerans",
    "Cyclopoid Copepods",
    "Larval Odonates"
  )

perch_relabund_long$treatment <-
  as.factor(perch_relabund_long$corral)

levels(perch_relabund_long$treatment) <-
  c("0(1)", "414", "29,240", "100", "6", "7,071", "0(2)", "1,710")

png(
  "Perch Diet Plot by Taxa Relative Abundance.png",
  width = 12,
  height = 8,
  units = "cm",
  res = 500
)

ggplot(subset(perch_relabund_long, total.animals != 0)) +
  geom_col(aes(
    x = ID,
    y = count,
    fill = reorder(taxa, 1 / (count + 1), mean)
  ),
  colour = "black",
  size = 0.25) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(paste("Dose (MPs" ~ L ^ -1 * ")")),
       y = "Proportion of Individuals") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Taxa") +
  facet_grid(
    . ~ reorder(treatment, MPconcentration, mean),
    scales = "free_x",
    switch = "x",
    space = "free_x"
  ) +
  theme1 +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm")
  )

dev.off()


# GLMM for total individuals ----

# Poisson
mod1 <- glmmTMB(
  total.animals ~
    log(MPconcentration + 6) +
    body.weight +
    YP.end,
  family = poisson(link = "log"),
  data = perch_diet
)
summary(mod1)
res1 <- simulateResiduals(mod1)
plot(res1)  # BAD

# NB with linear variance
mod2 <- glmmTMB(
  total.animals ~
    log(MPconcentration + 6) +
    body.weight +
    YP.end,
  family = nbinom1(link = "log"),
  data = perch_diet
)
summary(mod2)
res2 <- simulateResiduals(mod2)
plot(res2)  # DECENT

# NB with quadratic variance
mod3 <- glmmTMB(
  total.animals ~
    log(MPconcentration + 6) +
    body.weight +
    YP.end,
  family = nbinom2(link = "log"),
  data = perch_diet
)
summary(mod3)
res3 <- simulateResiduals(mod3)
plot(res3)  # BEST

## Plot model predictions ----

ingestionpred <-
  as.data.frame(ggemmeans(mod3,
                          terms = "YP.end")) %>%
  rename(YP.end = x)

png(
  "Perch Diet Totals.png",
  width = 9,
  height = 7,
  units = "cm",
  res = 600
)

ggplot() +
  geom_ribbon(
    data = ingestionpred,
    aes(x = YP.end,
        ymin = conf.low,
        ymax = conf.high),
    fill = "lime green",
    alpha = 0.3
  ) +
  geom_line(data = ingestionpred,
            aes(x = YP.end,
                y = predicted),
            size = 0.5) +
  geom_point(data = perch_diet,
             aes(x = YP.end,
                 y = total.animals)) +
  scale_y_continuous(expand = c(0.009, 0.005),
                     breaks = c(0, 1, 10, 100, 1000)) +
  coord_trans(y = "log1p") +
  labs(x = "Number of Surviving Perch",
       y = "Total number of Individuals (log scale)") +
  theme1

dev.off()

# Interestingly in mesocosm with more surviving perch, individuals seemed to
# be eating a higher number of animals


# Check that total individuals is a good measure of gut fullness ----

gutmod1 <-
  lm(log(gi.weight) ~ log(total.animals + 1), data = perch_diet)

plot(resid(gutmod1) ~ fitted(gutmod1), data = perch_diet)
abline(0, 0)

gutmodresid <- simulateResiduals(gutmod1)
plot(gutmodresid)  # looks good

summary(gutmod1)
# total indv. is a decent predictor of GI weight (R2 = 0.35)

gutpred <- predict(gutmod1, se.fit = TRUE)

gutpred$lower <- with(gutpred, fit - se.fit)
gutpred$upper <- with(gutpred, fit + se.fit)


## Plot model predictions ----

png(
  "Fullness Plot.png",
  width = 9,
  height = 8,
  units = "cm",
  res = 500
)

ggplot(data = perch_diet) +
  geom_point(aes(x = total.animals,
                 y = gi.weight)) +
  geom_ribbon(
    aes(
      x = total.animals,
      ymin = exp(gutpred$lower),
      ymax = exp(gutpred$upper)
    ),
    fill = "aquamarine",
    alpha = 0.3
  ) +
  geom_line(aes(x = total.animals,
                y = exp(gutpred$fit)),
            size = 0.5,
            linetype = "dashed") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  labs(x = "Total Number of Individuals (log scale)",
       y = "Gastrointestinal Tract Weight (g)") +
  theme1

dev.off()

# Multivariate comparison ----

## Prepare data ----

# Response matrix

Y <- perch_diet[12:17]

# Predictor matrix

X <- perch_diet[, -c(12:17)]

## NMDS ----

# Remove rows with all zeros
sum <- numeric()
for (i in 1:nrow(Y)) {
  sum[i] <- sum(Y[i, ])
}

Y2 <- Y[sum > 0, ]
X2 <- X[sum > 0, ]

set.seed(2454)

dietnMDS <- metaMDS(Y2,
                    distance = "bray")

dietnMDS

dietnMDS$stress

# Shepards test/goodness of fit
stressplot(dietnMDS)

gof <- goodness(dietnMDS)
plot(dietnMDS, type = "t")
points(dietnMDS, display = "sites", cex = gof * 300)

# Site scores

diet_nmds_site_scores <- as.data.frame(scores(dietnMDS,
                                              display = "site"))
diet_nmds_site_scores <- cbind(diet_nmds_site_scores,
                               X2)

diet_nmds_site_scores$label <-
  diet_nmds_site_scores$corral

levels(diet_nmds_site_scores$label) <-
  c(c("0(1)",
      "414",
      "29,240",
      "100",
      "6",
      "7,071",
      "0(2)",
      "1,710"))

# Species scores

diet_nmds_species_scores <- as.data.frame(scores(dietnMDS,
                                                 display = "species"))

diet_nmds_species_scores$Taxa <- c(
  "Chironomid Pupae",
  "Amphipods",
  "Larval Chironomids",
  "Cladocerans",
  "Cyclopoid Copepods",
  "Larval Odonates"
)

### Generate hulls ----

diet_nmds_hulls <- data.frame()

for (i in 1:length(unique(diet_nmds_site_scores$corral))) {
  hull <-
    diet_nmds_site_scores[diet_nmds_site_scores$corral ==
                            unique(diet_nmds_site_scores$corral)[i],][chull(diet_nmds_site_scores[diet_nmds_site_scores$corral ==
                                                                                                    unique(diet_nmds_site_scores$corral)[i],
                                                                                                  c(1:2)]),]
  diet_nmds_hulls <- rbind(diet_nmds_hulls, hull)
}

diet_nmds_hulls$label <-
  diet_nmds_hulls$corral

levels(diet_nmds_hulls$label) <-
  c(c("0(1)",
      "414",
      "29,240",
      "100",
      "6",
      "7,071",
      "0(2)",
      "1,710"))

### Plot ----

png(
  "Perch nMDS.png",
  width = 8,
  height = 8,
  units = "cm",
  res = 500
)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_polygon(
    data = diet_nmds_hulls,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = reorder(label,
                     MPconcentration,
                     mean),
      colour = reorder(label,
                       MPconcentration,
                       mean)
    ),
    alpha = 0.3,
    size = 0.5
  ) +
  geom_point(
    data = diet_nmds_site_scores,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = reorder(label,
                     MPconcentration,
                     mean)
    ),
    size = 2,
    alpha = 0.75,
    shape = 21
  ) +
  geom_segment(
    data = diet_nmds_species_scores,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    alpha = 0.5,
    linewidth = 0.75,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "blue3"
  ) +
  geom_text(
    data = diet_nmds_species_scores,
    aes(
      x = 1.2 * NMDS1,
      y = 1.2 * NMDS2,
      label = Taxa
    ),
    size = 7 / .pt,
    colour =  "blue3"
  ) +
  scale_fill_viridis_d(name = "Corral",
                       option = "inferno",
                       direction = -1) +
  scale_colour_viridis_d(name = "Corral",
                         option = "inferno",
                         direction = -1) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 1)) +
  theme1 +
  theme(legend.position = "bottom")

dev.off()

## PERMANOVA ----

set.seed(425)

diet_PERMANOVA <-
  adonis2(
    Y2 ~ log(MPconcentration + 6) + body.weight,
    strata = X2$corral,
    method = "bray",
    by = "margin",
    data = X2
  )

diet_PERMANOVA

# Summarize diet by corral

diet_summary <-
  perch_diet_long %>%
  group_by(corral, taxa) %>%
  summarize(mean = mean(count)) %>%
  ungroup() %>%
  pivot_wider(values_from = mean,
              names_from = taxa)

# Look at biometric effects on Copepod ingestion ----

perch_diet$stdate <-
  as.numeric(as.Date(perch_diet$collection.date,
                     format = "%d-%m-%y")) -
  min(as.numeric(as.Date(perch_diet$collection.date,
                         format = "%d-%m-%y")))

copepodmod1 <-
  glmmTMB(
    cyclopoida ~
      body.length +
      log(MPconcentration + 6) +
      (1 | corral),
    family = nbinom1(link = "log"),
    data = perch_diet
  )

plot(simulateResiduals(copepodmod1))

plotResiduals(copepodmod1, perch_diet$stdate)

summary(copepodmod1)

## Basically can't get a well-fitting model from the data
