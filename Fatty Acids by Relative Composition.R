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

# Load data ----

FAs_percent <- read.csv("FAs_percent.csv", header = TRUE)
str(FAs_percent)
FAs_percent$ID <- as.factor(FAs_percent$ID)
FAs_percent$corral <- as.factor(FAs_percent$corral)
FAs_percent$sample.type <- as.factor(FAs_percent$sample.type)
FAs_percent$date <- as.Date(FAs_percent$date,
                            format = "%b. %d, %Y")

treatments <- data.frame(corral =
                           as.factor(c("B", "C", "D", "E", "F", "G", "H", "I")),
                         MPconcentration = as.numeric(c(0, 414, 29240, 100, 6,
                                                        7071, 0, 1710)))

FAs_percent <- left_join(FAs_percent,
                         treatments,
                         by = "corral")

# Prepare data ----

# Assign corral A MP concentration of 24 for zoop_FA

FAs_percent$MPconcentration[is.na(FAs_percent$MPconcentration)] <-
  24

# Convert to proportions
FAs_prop <- FAs_percent
FAs_prop[7:50] <- FAs_percent[7:50] / 100

## Separate out by sample type ----

perch_FA_prop <- subset(FAs_prop, sample.type == "perch")
zoop_FA_prop <- subset(FAs_prop, sample.type == "zooplankton")
food_FA_prop <- subset(FAs_prop, sample.type == "fish food")

## Prep perch data ----

### Add in perch biometrics data ----

perch_biometrics <- read.csv("perch_biometrics.csv", header = TRUE)

str(perch_biometrics)

perch_biometrics$corral <- as.factor(perch_biometrics$corral)
perch_biometrics$ID <- as.factor(perch_biometrics$ID)
perch_biometrics$sex <- as.factor(perch_biometrics$sex)

# Combine with perch FA data

perch_FA_prop2 <- left_join(perch_FA_prop,
                            perch_biometrics,
                            by = c("ID", "corral"))

### Add in perch population data ----

perch_pop <- read.csv("fish_pop.csv", header = TRUE)

# Combine with perch FA data

perch_FA_prop2 <- left_join(perch_FA_prop2,
                            perch_pop,
                            by = c("corral", "MPconcentration"))

## Prep zooplankton data ----

# Sample comparison ----

## Put data into long form ----

FA_prop_long <-
  FAs_prop[, c(1:3, 5, 6, 17, 20, 23, 26, 27, 31, 35, 36, 40, 42, 45:47, 51)] %>%
  pivot_longer(names(FAs_prop[c(17, 20, 23, 26, 27, 31, 35, 36, 40, 42, 45:47)]),
               names_to = "metric",
               values_to = "value")

FA_prop_long$metric <- as.factor(FA_prop_long$metric)

levels(FA_prop_long$sample.type) <-
  c("Mysis Fish Food",
    "Yellow Perch",
    "Zooplankton")

# Focus on endpoint Zooplankton
FA_prop_long <- subset(FA_prop_long,
                       date >= "2021-08-01" |
                         sample.type == "Mysis Fish Food")

# Set MP concentration to 0 for fish food

FA_prop_long$MPconcentration[FA_prop_long$sample.type ==
                               "Mysis Fish Food"] <- 0

# Define labeller

prop_labels <- as_labeller(
  c(
    "C_16.1n.7" = "Palmitoleic acid",
    "C_18.1n.9" = "Oleic acid",
    "C_18.2n.6" = "Linoleic acid",
    "C_18.3n.3" = "Alpha-linoleic acid",
    "C_20.4n.6" = "Arachidonic acid",
    "C_20.5n.3" = "Eicosapentaenoic acid",
    "C_22.6n.3" = "Docosahexaenoic acid",
    "HUFAs" = "Total HUFAs",
    "PUFAs" = "Total PUFAs",
    "total_MUFAs" = "Total MUFAs",
    "total_N.3_PUFAs" = "Total n-3 PUFAs",
    "total_N.6_PUFAs" = "Total n-6 PUFAs",
    "total_SFAs" = "Total SFAs"
  )
)

## Plot ----

png(
  "Sample Comparison FA Plot.png",
  width = 19,
  height = 19,
  units = "cm",
  res = 600
)

ggplot(FA_prop_long) +
  geom_point(aes(x = MPconcentration,
                 y = value,
                 colour = sample.type)) +
  geom_smooth(aes(x = MPconcentration,
                  y = value,
                  colour = sample.type),
              method = "lm") +
  facet_wrap( ~ metric,
              scales = "free_y",
              labeller = prop_labels,
              ncol = 3) +
  labs(x = "Sample Type",
       y = "Proportion of Total Fatty Acids") +
  scale_colour_manual(values = c("orange", "purple", "blue"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

dev.off()


# Analyses ----

## Perch analyses ----

### Individual FAs ----

#### DHA ----

DHA.prop.mod1 <-
  glmmTMB(C_22.6n.3 ~
            body.weight +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA_prop2)

plot(simulateResiduals(DHA.prop.mod1))

summary(DHA.prop.mod1)  # no effect of MP concentration

#### EPA ----

EPA.prop.mod1 <-
  glmmTMB(C_20.5n.3 ~
            body.weight +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA_prop2)

plot(simulateResiduals(EPA.prop.mod1))

summary(EPA.prop.mod1)  # positive correlation with body size

# Plot model predictions

EPA.prop.pred <-
  ggemmeans(EPA.prop.mod1,
            terms = "body.weight") %>%
  rename(body.weight = x)

png(
  "Perch Compositional EPA Plot.png",
  width = 8.84,
  height = 5,
  units = "cm",
  res = 500
)

ggplot() +
  geom_ribbon(
    data = EPA.prop.pred,
    aes(x = body.weight,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.3
  ) +
  geom_line(data = EPA.prop.pred,
            aes(x = body.weight,
                y = predicted)) +
  geom_point(data = perch_FA_prop2,
             aes(x = body.weight,
                 y = C_20.5n.3),
             size = 1) +
  labs(x = "Body Weight (g)",
       y = "Proportion EPA") +
  theme1

dev.off()

#### ARA ----

ARA.prop.mod1 <-
  glmmTMB(C_20.4n.6 ~
            body.weight +
            log(MPconcentration + 6) +
            (1 | corral),
          data = perch_FA_prop2)

plot(simulateResiduals(ARA.prop.mod1))

summary(ARA.prop.mod1)  # no effect of MP concentration


### Multivariate analyses ----

#### Pull out FA data and covariates ----

# Remove any FAs that contribute less than 1% and re-scale response so it
# sums to 1

trimmed_perch_FA <-
  data.frame(perch_FA_prop2[, c(7:16, 18:25, 27:34, 36:44)] %>%
               select(where(function(x) {
                 mean(x) >= 0.01
               })))

perch_FA_prop_covariates <- perch_FA_prop2[, c(1:3, 51:59)]

perch_FA_prop_totals <- perch_FA_prop2[, c(17, 26, 35, 45)]

summary(perch_FA_prop_totals)

trimmed.FA.names <-
  c(
    "14:0",
    "16:0",
    "18:0",
    "16:1n-7",
    "18:1n-7",
    "18:1n-9",
    "18:2n-6",
    "20:4n-6",
    "22:5n-6",
    "18:3n-3",
    "18:4n-3",
    "20:5n-3",
    "22:5n-3",
    "22:6n-3"
  )

#### CCA ----
set.seed(424)

perch_FA_prop_cca <-
  cca((trimmed_perch_FA) ~
        log(MPconcentration + 6) +
        body.weight,
      data = perch_FA_prop_covariates)

summary(perch_FA_prop_cca)

RsquareAdj(perch_FA_prop_cca)

screeplot(perch_FA_prop_cca)

anova(perch_FA_prop_cca, by = "term")
anova(perch_FA_prop_cca, by = "margin")

plot(perch_FA_prop_cca,
     scaling = 1,
     display = c("lc", "wa", "cn"))
plot(perch_FA_prop_cca,
     scaling = 2,
     display = c("sp", "cn"))


#### Pull out scores for scaling 1----

# 'Site' scores
perch_FA_prop_cca_site1 <-
  as.data.frame(scores(perch_FA_prop_cca, display = "lc",
                       scaling = 1))

perch_FA_prop_cca_site1 <- cbind(perch_FA_prop_covariates,
                                 perch_FA_prop_cca_site1[, 1:2])
perch_FA_prop_cca_bp1 <-
  as.data.frame(scores(perch_FA_prop_cca,
                       display = "bp",
                       scaling = 1))

perch_FA_prop_cca_bp1$label <- c("MP Concentration", "Body Weight")

#### Pull out scores for scaling 2----

# 'Species' scores
perch_FA_prop_cca_species2 <-
  as.data.frame(scores(perch_FA_prop_cca, display = "species",
                       scaling = 2))

perch_FA_prop_cca_species2$FA <- trimmed.FA.names

perch_FA_prop_cca_bp2 <-
  as.data.frame(scores(perch_FA_prop_cca,
                       display = "bp",
                       scaling = 2))

perch_FA_prop_cca_bp2$label <- c("MP Concentration", "Body Weight")

#### Plot scaling 1 ----

perchccaplot1 <-
  ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_point(
    data = perch_FA_prop_cca_site1,
    aes(
      x = CCA1,
      y = CCA2,
      fill = as.factor(MPconcentration)
    ),
    size = 2,
    shape = 21,
    alpha = 0.95,
    colour = "grey30"
  ) +
  geom_text(
    data = perch_FA_prop_cca_bp1,
    aes(
      x =  CCA1 * 1.1,
      y = CCA2 * 1.1,
      label = label
    ),
    size = 10 / .pt,
    colour =  "grey10"
  ) +
  geom_segment(
    data = perch_FA_prop_cca_bp1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2,
      colour = vars
    ),
    alpha = 0.75,
    linewidth = 1,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "grey10"
  ) +
  scale_colour_viridis_d(name =
                           expression(paste("Exposure Concentration (MPs" ~
                                              L ^ -1 * ")")),
                         option = "plasma") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  scale_x_continuous(limits = c(-0.35, 0.35)) +
  labs(x = "CCA1",
       y = "CCA2") +
  theme1 +
  theme(
    legend.key.size = unit(0.2, "cm"),
    legend.spacing = unit(0, "cm"),
    legend.position = "bottom"
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.position = "top"
    ),
    colour = guide_legend(
      nrow = 2,
      byrow =  TRUE,
      title.position = "top"
    )
  )


#### Plot scaling 2 ----

perchccaplot2 <-
  ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_text_repel(
    data = perch_FA_prop_cca_species2,
    aes(x = CCA1,
        y = CCA2,
        label = FA),
    size = 8 / .pt,
    colour = "blue3",
    box.padding = 0
  ) +
  geom_text(
    data = perch_FA_prop_cca_bp2,
    aes(
      x = CCA1 * 0.22,
      y = CCA2 * 1.25,
      label = label
    ),
    size = 10 / .pt,
    colour =  "grey10"
  ) +
  geom_segment(
    data = perch_FA_prop_cca_bp2,
    aes(
      x = 0,
      y = 0,
      xend = CCA1 * 0.3,
      yend = CCA2 * 0.3,
      colour = vars
    ),
    alpha = 0.75,
    linewidth = 0.5,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "grey10"
  ) +
  scale_colour_viridis_d(name =
                           expression(paste("Exposure Concentration (MPs" ~
                                              L ^ -1 * ")")),
                         option = "plasma") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  labs(x = "CCA1",
       y = "CCA2") +
  theme1 +
  theme(legend.position = "none")

#### Combine plots ----

png(
  "Perch FA Proportions CCA.png",
  width = 18,
  height = 10,
  units = "cm",
  res = 600
)

plot_grid(
  perchccaplot1,
  perchccaplot2,
  ncol = 2,
  align = "h",
  axis = "b",
  labels = c("A", "B"),
  label_size = 12
)

dev.off()


#### Dive in factors driving variation ----

trimmed_perch_FA3 <-
  trimmed_perch_FA[!is.na(perch_FA_prop_covariates$sex) &
                     !is.na(perch_FA_prop_covariates$gonad.weight), ]
perch_FA_prop_covariates2 <-
  perch_FA_prop_covariates[!is.na(perch_FA_prop_covariates$sex) &
                             !is.na(perch_FA_prop_covariates$gonad.weight), ]

perch_FA_prop_covariates2$corral <-
  as.factor(perch_FA_prop_covariates2$corral)

set.seed(4554)

perch_FA_prop_cca2 <-
  cca(trimmed_perch_FA3 ~
        log(MPconcentration + 6) +
        body.weight +
        gonad.weight +
        sex,
      data = perch_FA_prop_covariates2)

summary(perch_FA_prop_cca2)

RsquareAdj(perch_FA_prop_cca2)

screeplot(perch_FA_prop_cca2)

anova(perch_FA_prop_cca2, by = "term")  # gonad weight sig
anova(perch_FA_prop_cca2, by = "margin")

plot(perch_FA_prop_cca2,
     scaling = 1,
     display = c("wa", "cn"))
plot(perch_FA_prop_cca2,
     scaling = 2,
     display = c("sp", "cn"))

###### Pull out scores for scaling 2----

# 'Site' scores
perch_FA_prop_cca_site3 <-
  as.data.frame(scores(perch_FA_prop_cca2, display = "site",
                       scaling = 2))

# 'Species' scores
perch_FA_prop_cca_species3 <-
  as.data.frame(scores(perch_FA_prop_cca2, display = "species",
                       scaling = 2))

perch_FA_prop_cca_species3$FA <- trimmed.FA.names

perch_FA_prop_cca_site3 <- cbind(perch_FA_prop_covariates2,
                                 perch_FA_prop_cca_site3[, 1:2])

perch_FA_prop_cca_centroids3 <-
  as.data.frame(scores(perch_FA_prop_cca2, display = "cn",
                       scaling = 2))

perch_FA_prop_cca_bp3 <-
  as.data.frame(scores(perch_FA_prop_cca2,
                       display = "bp",
                       scaling = 2))

# select just gonad weight
perch_FA_prop_cca_bp3 <- perch_FA_prop_cca_bp3[3, ]

perch_FA_prop_cca_bp3$label <- "Gonad Weight"


##### Plot ----

png(
  "Perch FA Proportions CCA Reduced.png",
  width = 8.84,
  height = 7.8,
  units = "cm",
  res = 600
)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_text_repel(
    data = perch_FA_prop_cca_species3,
    aes(x = CCA1,
        y = CCA2,
        label = FA),
    size = 8 / .pt,
    colour = "blue3",
    box.padding = 0,
    max.overlaps = 20,
    alpha = 0.75
  ) +
  geom_text(
    data = perch_FA_prop_cca_bp3,
    aes(
      x = CCA1 * 1.1,
      y = CCA2 * 1.1,
      label = label
    ),
    size = 8 / .pt,
    colour =  "grey10"
  ) +
  geom_segment(
    data = perch_FA_prop_cca_bp3,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2,
      colour = vars
    ),
    alpha = 0.75,
    linewidth = 0.5,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "grey10"
  ) +
  scale_colour_viridis_d(name =
                           expression(paste("Exposure Concentration (MPs" ~
                                              L ^ -1 * ")")),
                         option = "plasma") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  labs(x = "CCA1",
       y = "CCA2") +
  scale_x_continuous(limits = c(-0.5, 1.1)) +
  theme1 +
  theme(legend.position = "none")

dev.off()

# These are the plots with everything

plot(perch_FA_prop_cca,
     scaling = 1,
     display = c("lc", "cn"))

plot(perch_FA_prop_cca,
     scaling = 2,
     display = c("sp", "cn"))

summary(perch_FA_prop_cca)


### Explore individual patterns ----

#### Pairwise correlations ----

pairs(trimmed_perch_FA)

#### Look at short- vs. long-chain PUFAs ----

trimmed_perch_FA2 <-
  cbind(trimmed_perch_FA, perch_FA_prop_covariates)

trimmed_perch_FA2_long <-
  pivot_longer(trimmed_perch_FA2,
               cols = c(1:14),
               names_to = "FA")

trimmed_perch_FA2_long$FA <- as.factor(trimmed_perch_FA2_long$FA)

levels(trimmed_perch_FA2_long$FA)

reduced_perch_FA2_long <-
  trimmed_perch_FA2_long %>%
  filter(
    FA == "C_18.1n.9" |
      FA == "C_18.2n.6" |
      FA == "C_18.3n.3" |
      FA == "C_18.4n.3" |
      FA == "C_20.4n.6" |
      FA == "C_22.5n.6" |
      FA == "C_22.6n.3"
  )

ggplot(reduced_perch_FA2_long) +
  geom_point(aes(x = body.weight,
                 y = value,
                 colour = FA)) +
  geom_smooth(aes(x = body.weight,
                  y = value,
                  colour = FA),
              method = "lm") +
  facet_grid(sex ~ FA, scale = "free") +
  theme1

ggplot(reduced_perch_FA2_long) +
  geom_point(aes(x = gonad.weight,
                 y = value,
                 colour = FA)) +
  geom_smooth(aes(x = gonad.weight,
                  y = value,
                  colour = FA),
              method = "lm") +
  facet_grid(sex ~ FA, scale = "free") +
  theme1

ggplot(reduced_perch_FA2_long) +
  geom_point(aes(x = liver.weight,
                 y = value,
                 colour = FA)) +
  geom_smooth(aes(x = liver.weight,
                  y = value,
                  colour = FA),
              method = "lm") +
  facet_wrap( ~ FA, scale = "free") +
  theme1

# Hard to figure out what's going on here


## Zooplankton analyses ----

### Individual FAs ----

# Separate just midpoint and endpoint zoops

zoop.end <-
  zoop_FA_prop %>%
  filter(date == "2021-08-09" |
           date == "2021-07-06")

#### DHA ----

DHA.prop.zoop.mod1 <-
  glmmTMB(C_22.6n.3 ~
            log(MPconcentration + 6) +
            as.factor(date) +
            (1 | corral),
          data = zoop.end)

plot(simulateResiduals(DHA.prop.zoop.mod1))

summary(DHA.prop.zoop.mod1)
# no effect of MP concentration
# more DHA at the end

#### EPA ----

EPA.prop.zoop.mod1 <-
  glmmTMB(C_20.5n.3 ~
            log(MPconcentration + 6) +
            as.factor(date) +
            (1 | corral),
          data = zoop.end)

plot(simulateResiduals(EPA.prop.zoop.mod1))

summary(EPA.prop.zoop.mod1)
# no effect of MP concentration
# less EPA at the end

#### ARA ----

ARA.prop.zoop.mod1 <-
  glmmTMB(C_20.4n.6 ~
            MPconcentration +
            as.factor(date) +
            (1 | corral),
          data = zoop.end)

plot(simulateResiduals(ARA.prop.zoop.mod1))

summary(ARA.prop.zoop.mod1)
# no effect of MP concentration
# less ARA at the end


#### Plot HUFA composition ----

zoop.end.long <-
  zoop.end %>%
  select(-c(7:30, 32:39, 41, 43:50)) %>%
  pivot_longer(cols = c(7:9),
               names_to = "FA")

zoop.end.long$timepoint <-
  as.factor(zoop.end.long$date) %>%
  recode_factor(`2021-07-06` = "Day 34",
                `2021-08-09` = "Day 68")

png(
  "Zoop Compositional HUFA Plot.png",
  width = 8.84,
  height = 6,
  units = "cm",
  res = 500
)

ggplot(zoop.end.long) +
  geom_boxplot(aes(x = timepoint,
                   y = value,
                   fill = FA),
               alpha = 0.6) +
  scale_fill_viridis_d(
    name = "",
    labels = c("ARA",
               "EPA",
               "DHA"),
    option = "magma"
  ) +
  labs(x = "Date",
       y = "Proportion") +
  theme1

dev.off()

### Multivariate ----

#### Add zooplankton biomass data ----

zoopbm <- read.csv("zoop_biomass_2021.csv",
                   header = TRUE,
                   stringsAsFactors = TRUE)

str(zoopbm)

# Match the time points

unique(zoop_FA_prop$date)

zoop_FA_prop <-
  zoop_FA_prop %>%
  mutate(timepoint = as.factor(as.character(date)))

levels(zoop_FA_prop$timepoint) <- c("Beginning",
                                    "Mid-point",
                                    "End")

unique(zoopbm$date_sampled)

zoopbm2 <-
  zoopbm %>%
  filter(
    date_sampled == "2021-06-01" |
      date_sampled == "2021-07-05" |
      date_sampled == "2021-08-09"
  ) %>%
  rename(timepoint = date_sampled)

zoopbm2$timepoint <- as.character(zoopbm2$timepoint)
zoopbm2$timepoint <- as.factor(zoopbm2$timepoint)

levels(zoopbm2$timepoint) <- c("Beginning", "Mid-point", "End")

zoopbm2 <- zoopbm2[, -1]

zoopbm2 <-
  zoopbm2 %>%
  pivot_wider(names_from = "order",
              values_from = "biomass_ug_l")

zoop_FA_prop2 <-
  left_join(zoop_FA_prop,
            zoopbm2,
            by = c("corral",
                   "timepoint"))

#### Pull out covariates ----

# Remove any FAs that contribute less than 1% and re-scale response so it
# sums to 1

trimmed_zoop_FA <-
  data.frame(zoop_FA_prop2[, c(7:16, 18:25, 27:34, 36:44)]) %>%
  select(where(function(x) {
    mean(x) >= 0.01
  }))

zoop_FA_prop_covariates <- zoop_FA_prop2[, c(1:3, 51:56), ]

zoop_FA_prop_totals <- zoop_FA_prop2[, c(17, 26, 35, 45)]

summary(zoop_FA_prop_totals)

trimmed.FA.names.zoops <-
  c(
    "14:0",
    "16:0",
    "18:0",
    "16:1n-7",
    "18:1n-7",
    "18:1n-9",
    "22:1n-9",
    "18:2n-6",
    "18:3n-6",
    "20:4n-6",
    "22:5n-6",
    "18:3n-3",
    "18:4n-3",
    "20:5n-3",
    "22:6n-3"
  )

#### CCA ----

set.seed(4399)

zoop_FA_prop_cca <-
  cca(
    trimmed_zoop_FA ~
      log(MPconcentration + 6) +
      as.factor(date) +
      Calanoida +
      Cladocera +
      Cyclopoida,
    scale. = FALSE,
    data = zoop_FA_prop_covariates
  )

RsquareAdj(zoop_FA_prop_cca)

plot(zoop_FA_prop_cca,
     scaling = 1,
     display = c("lc", "cn"))

plot(zoop_FA_prop_cca,
     scaling = 2,
     display = c("sp", "cn"))

screeplot(zoop_FA_prop_cca)

summary(zoop_FA_prop_cca)

anova(zoop_FA_prop_cca, by = "term")
anova(zoop_FA_prop_cca, by = "margin")
anova(zoop_FA_prop_cca, by = "onedf")

#### Pull out scaling 1 scores ----

# 'Site' scores
zoop_FA_prop_cca_site1 <-
  as.data.frame(scores(zoop_FA_prop_cca, display = "site",
                       scaling = 1))

# 'Species' scores
zoop_FA_prop_cca_species1 <-
  as.data.frame(scores(zoop_FA_prop_cca, display = "species",
                       scaling = 1))

zoop_FA_prop_cca_species1$FA <- trimmed.FA.names.zoops

zoop_FA_prop_cca_site1 <- cbind(zoop_FA_prop_covariates,
                                zoop_FA_prop_cca_site1[, 1:2])

zoop_FA_prop_cca_site1$label <- zoop_FA_prop_cca_site1$corral

levels(zoop_FA_prop_cca_site1$label) <- c("24",
                                          "0(1)",
                                          "414",
                                          "29,240",
                                          "100",
                                          "6",
                                          "7,071",
                                          "0(2)",
                                          "1,710")

zoop_FA_prop_cca_site1$label <-
  factor(
    zoop_FA_prop_cca_site1$label,
    levels = c(
      "0(1)",
      "0(2)",
      "6",
      "24",
      "100",
      "414",
      "1,710",
      "7,071",
      "29,240"
    )
  )


zoop_FA_prop_cca_site1$point <-
  as.factor(as.character(zoop_FA_prop_cca_site1$date))

levels(zoop_FA_prop_cca_site1$point) <- c("Day -6",
                                          "Day 34",
                                          "Day 68")
# Generate hulls
zoop_FA_date_hulls <- data.frame()

for (i in 1:length(unique(zoop_FA_prop_cca_site1$point))) {
  subset <-
    zoop_FA_prop_cca_site1[zoop_FA_prop_cca_site1$point ==
                             unique(zoop_FA_prop_cca_site1$point)[i],]
  hullrows <- chull(subset[, c(8:9)])
  hull <- subset[hullrows, ]
  zoop_FA_date_hulls <- rbind(zoop_FA_date_hulls, hull)
}

zoop_FA_prop_cca_cn1 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 1))

zoop_FA_prop_cca_cn1$label <-
  factor(c("Day -6",
           "Day 34",
           "Day 68"))

zoop_FA_prop_cca_bp1 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "bp",
                    scaling = 1))[c(4:6), ]

zoop_FA_prop_cca_bp1$label = c("Calanoida",
                               "Cladocera",
                               "Cyclopoida")

##### Plot ----

zoopFAplot1 <-
  ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_polygon(
    data = zoop_FA_date_hulls,
    aes(x = CCA1,
        y = CCA2,
        group = point),
    alpha = 0.25,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    data = zoop_FA_prop_cca_site1,
    aes(
      x = CCA1,
      y = CCA2,
      fill = label,
      shape = point
    ),
    size = 3,
    alpha = 0.75,
    colour = "grey30"
  ) +
  geom_text(
    data = zoop_FA_prop_cca_bp1,
    aes(
      x =  CCA1 * 1.4,
      y = CCA2 * 1.4,
      label = label
    ),
    size = 8 / .pt,
    colour =  "grey10"
  ) +
  geom_segment(
    data = zoop_FA_prop_cca_bp1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2,
      colour = vars
    ),
    alpha = 0.5,
    linewidth = 0.75,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "grey10"
  ) +
  geom_text_repel(
    data = zoop_FA_prop_cca_cn1,
    aes(x = CCA1,
        y = CCA2,
        label = label),
    size = 10 / .pt,
    colour =  "blue3",
    box.padding = unit(0, "cm")
  ) +
  scale_shape_manual(values = c(21, 24, 22),
                     name = "Time Point") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "inferno",
                       direction = -1) +
  labs(x = "CCA1",
       y = "CCA2") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme1 +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0, "cm"))

#### Pull out scaling 2 scores ----

# 'Site' scores
zoop_FA_prop_cca_site2 <-
  as.data.frame(scores(zoop_FA_prop_cca, display = "site",
                       scaling = 2))

# 'Species' scores
zoop_FA_prop_cca_species2 <-
  as.data.frame(scores(zoop_FA_prop_cca, display = "species",
                       scaling = 2))

zoop_FA_prop_cca_species2$FA <- trimmed.FA.names.zoops

zoop_FA_prop_cca_site2 <- cbind(zoop_FA_prop_covariates,
                                zoop_FA_prop_cca_site1[, 1:2])

zoop_FA_prop_cca_site2$label <- zoop_FA_prop_cca_site2$corral

levels(zoop_FA_prop_cca_site2$label) <- c("24",
                                          "0(1)",
                                          "414",
                                          "29,240",
                                          "100",
                                          "6",
                                          "7,071",
                                          "0(2)",
                                          "1,710")

zoop_FA_prop_cca_site2$label <-
  factor(
    zoop_FA_prop_cca_site2$label,
    levels = c(
      "0(1)",
      "0(2)",
      "6",
      "24",
      "100",
      "414",
      "1,710",
      "7,071",
      "29,240"
    )
  )


zoop_FA_prop_cca_site2$point <-
  as.factor(as.character(zoop_FA_prop_cca_site2$date))

levels(zoop_FA_prop_cca_site2$point) <- c("Day -6",
                                          "Day 34",
                                          "Day 68")

zoop_FA_prop_cca_cn1 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 1))

zoop_FA_prop_cca_cn1$label <-
  factor(c("Day -6",
           "Day 34",
           "Day 68"))

zoop_FA_prop_cca_cn2 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 2))

zoop_FA_prop_cca_cn2$label <-
  factor(c("Day -6",
           "Day 34",
           "Day 68"))

zoop_FA_prop_cca_bp2 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "bp",
                    scaling = 2))[c(4:6), ]

zoop_FA_prop_cca_bp2$label = c("Calanoida",
                               "Cladocera",
                               "Cyclopoida")

##### Plot ----

zoopFAplot2 <-
  ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_text_repel(
    data = zoop_FA_prop_cca_species2,
    aes(x = CCA1,
        y = CCA2,
        label = FA),
    size = 8 / .pt,
    colour = "blue3",
    box.padding = 0,
    max.overlaps = 20
  ) +
  geom_text(
    data = zoop_FA_prop_cca_bp2,
    aes(
      x =  CCA1 * 1.4,
      y = CCA2 * 1.4,
      label = label
    ),
    size = 8 / .pt,
    colour =  "grey10"
  ) +
  geom_segment(
    data = zoop_FA_prop_cca_bp2,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2,
      colour = vars
    ),
    alpha = 0.5,
    linewidth = 0.75,
    arrow = arrow(
      angle = 20,
      length = unit(0.25, "cm"),
      type = "open"
    ),
    colour = "grey10"
  ) +
  geom_text_repel(
    data = zoop_FA_prop_cca_cn2,
    aes(x = CCA1,
        y = CCA2,
        label = label),
    size = 10 / .pt,
    colour =  "purple3",
    box.padding = unit(0, "cm")
  ) +
  scale_y_continuous(limits = c(-1.4, 1.1)) +
  scale_x_continuous(limits = c(-1.5, 1.1)) +
  scale_shape_manual(values = c(21, 24, 22),
                     name = "Time Point") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  labs(x = "CCA1",
       y = "CCA2") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme1 +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0, "cm"))

#### Combine plot ----

png(
  "Zooplankton FA Proportions CCA.png",
  width = 18,
  height = 15,
  units = "cm",
  res = 600
)

plot_grid(
  zoopFAplot1,
  zoopFAplot2,
  rel_widths = c(1, 0.6),
  ncol = 1,
  align = "v",
  axis = "r",
  labels = c("A", "B"),
  label_size = 12
)

dev.off()


# Plot proportions of major FA groups ----
## Fish food ----

food_FA_prop_long <-
  food_FA_prop %>%
  pivot_longer(cols = c(17, 26, 35, 45),
               names_to = "FA")

food_FA_prop_long$FA <- as.factor(food_FA_prop_long$FA)

levels(food_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

png(
  "Fish Food FAs Proportionss Plot.png",
  width = 8.84,
  height = 4,
  units = "cm",
  res = 600
)

ggplot(food_FA_prop_long) +
  geom_col(aes(x = ID,
               y = value,
               fill = FA)) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  labs(x = "Sample",
       y = "Proportion") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

dev.off()

## Perch ----

perch_FA_prop_long <-
  perch_FA_prop2 %>%
  pivot_longer(cols = c(17, 26, 35, 45),
               names_to = "FA")

perch_FA_prop_long$FA <- as.factor(perch_FA_prop_long$FA)

levels(perch_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

png(
  "Perch FAs Proportionss Plot.png",
  width = 18,
  height = 8.84,
  units = "cm",
  res = 600
)

ggplot(perch_FA_prop_long) +
  geom_col(aes(x = ID,
               y = value,
               fill = FA)) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  labs(x = "Sample",
       y = "Proportion") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  facet_wrap( ~ MPconcentration,
              scales = "free_x",
              nrow = 2) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

dev.off()

## Zooplankton ----

zoop_FA_prop_long <-
  zoop_FA_prop %>%
  pivot_longer(cols = c(17, 26, 35, 45),
               names_to = "FA")

zoop_FA_prop_long$FA <- as.factor(zoop_FA_prop_long$FA)

levels(zoop_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

zoop_FA_prop_long$timepoint <-
  as.factor(zoop_FA_prop_long$date)

levels(zoop_FA_prop_long$timepoint) <-
  c("Start",
    "Mid-point",
    "End")

png(
  "Zooplankton FAs Proportionss Plot.png",
  width = 18,
  height = 12,
  units = "cm",
  res = 600
)

ggplot(zoop_FA_prop_long) +
  geom_col(aes(
    x = reorder(ID, date),
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
       y = "Proportion") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  facet_wrap( ~ MPconcentration,
              scales = "free_x") +
  theme1 +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

dev.off()
