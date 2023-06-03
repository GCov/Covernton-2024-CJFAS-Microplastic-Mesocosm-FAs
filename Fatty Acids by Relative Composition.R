# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)
library(vegan)
library(tidyr)
library(DirichletReg)
library(easyCODA)
library(ggrepel)
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

FAs_percent$MPconcentration[is.na(FAs_percent$MPconcentration)] <- 24

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

### Center and scale variables ----

perch_FA_prop2$scaled.date <- 
  as.numeric(scale(as.numeric(perch_FA_prop2$date), 
                   center = min(as.numeric(perch_FA_prop2$date)),
                   scale = 9))

perch_FA_prop2$scaled.body.weight <-
  as.numeric(scale(perch_FA_prop2$body.weight,
                   center = TRUE))

perch_FA_prop2$scaled.MPconcentration <-
  as.numeric(scale(perch_FA_prop2$MPconcentration,
                   center = TRUE))

perch_FA_prop2$scaled.pop <-
  as.numeric(scale(perch_FA_prop2$YP.end,
                   center = TRUE))

# Re-scale response so it sums to 1

perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] 

adjusted_Y <- 
  DR_data(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] / 
            rowSums(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)]))

perch_FA_dir_mod1 <-
  DirichReg(adjusted_Y ~ scaled.body.weight + scaled.MPconcentration, 
            data = perch_FA_prop2)
perch_FA_dir_mod1
summary(perch_FA_dir_mod1)



## Prep zooplankton data ----

### Center and scale variables ----

zoop_FA_prop$scaled.date <- 
  as.numeric(scale(as.numeric(zoop_FA_prop$date), 
                   center = min(as.numeric(zoop_FA_prop$date)),
                   scale = 74))

zoop_FA_prop$scaled.MPconcentration <-
  as.numeric(scale(zoop_FA_prop$MPconcentration,
                   center = TRUE))

zoop_FA_prop$scaled.sample.mass.mg <-
  as.numeric(scale(zoop_FA_prop$sample.mass.mg,
                   center = TRUE))

# Sample comparison ----

## Put data into long form ----

FA_prop_long <- 
  FAs_prop[,c(1:3,5,6,17,20,23,26,27,31,35,36,40,42,45:47,51)] %>%
  pivot_longer(names(FAs_prop[c(17,20,23,26,27,31,35,36,40,42,45:47)]),
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

prop_labels <- as_labeller(c("C_16.1n.7" = "Palmitoleic acid",
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
                             "total_SFAs" = "Total SFAs"))

## Plot ----

png("Sample Comparison FA Plot.png",
    width = 19,
    height= 19, 
    units = "cm",
    res = 600)

ggplot(FA_prop_long) +
  geom_point(aes(x = MPconcentration,
                 y = value,
                 colour = sample.type)) +
  geom_smooth(aes(x = MPconcentration,
                  y = value,
                  colour = sample.type),
              method = "lm") +
  facet_wrap(~metric, scales = "free_y", labeller = prop_labels, ncol = 3) +
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
            TL + 
            log(MPconcentration + 1) +
            (1 | corral), 
          data = perch_FA_prop2)

plot(simulateResiduals(DHA.prop.mod1))

summary(DHA.prop.mod1)  # no effect of MP concentration

#### EPA ----

EPA.prop.mod1 <- 
  glmmTMB(C_20.5n.3 ~ 
            TL +
            log(MPconcentration + 1) +
            (1 | corral), 
          data = perch_FA_prop2)

  plot(simulateResiduals(EPA.prop.mod1))

summary(EPA.prop.mod1)  # positive correlation with body size

EPA.prop.pred <-
  ggemmeans(EPA.prop.mod1, 
            terms = "TL") %>% 
  rename(TL = x)

png("Perch Compositional EPA Plot.png",
    width = 8.84,
    height= 4, 
    units = "cm",
    res = 500)

ggplot() +
  geom_ribbon(data = EPA.prop.pred,
              aes(x = TL,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  geom_line(data = EPA.prop.pred,
            aes(x = TL,
                y = predicted)) +
  geom_point(data = perch_FA_prop2,
             aes(x = TL,
                 y = C_20.5n.3),
             size = 1) +
  labs(x = "Total Length (cm)",
       y = "Proportion EPA") +
  theme1

dev.off()
 
#### ARA ----

ARA.prop.mod1 <- 
  glmmTMB(C_20.4n.6 ~ 
            TL + 
            log(MPconcentration + 1) +
            (1 | corral), 
          data = perch_FA_prop2)

plot(simulateResiduals(ARA.prop.mod1))

summary(ARA.prop.mod1)  # no effect of MP concentration

### PCA ----

#### Pull out FA data and covariates ----

# Remove any FAs that contribute less than 1% and re-scale response so it 
# sums to 1

trimmed_perch_FA <- 
  data.frame(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] %>% 
               select(where(function(x){mean(x) >= 0.01})))

# trimmed_perch_FA <- 
#   trimmed_perch_FA / rowSums(trimmed_perch_FA)

perch_FA_prop_covariates <- perch_FA_prop2[,c(1:3,51:67),]

perch_FA_prop_totals <- perch_FA_prop2[,c(17,26,35,45)]

summary(perch_FA_prop_totals)

# Convert to CLRs

trimmed_perch_FA_clr <- decostand(trimmed_perch_FA,
                                  method = "clr",
                                  MARGIN = 1)

#### Run PCA ----
perch_FA_prop_pca <- rda(trimmed_perch_FA_clr,
                         scale. = FALSE)

# Bar plot of relative eigenvalues
barplot(as.vector(perch_FA_prop_pca$CA$eig)/sum(perch_FA_prop_pca$CA$eig))

# Calculate percentage of variance explained by first 2 aaxes
sum((as.vector(perch_FA_prop_pca$CA$eig)/sum(perch_FA_prop_pca$CA$eig))[1:2])

summary(perch_FA_prop_pca)

# 'Site' scores
perch_FA_prop_PCA_site <- perch_FA_prop_pca$CA$u

# 'Species' scores
perch_FA_prop_PCA_species <- data.frame(perch_FA_prop_pca$CA$v)

trimmed.FA.names <- 
  c("14:0", "16:0", "18:0", "16:1(n-7)", "18:(1n-7)", "18:1(n-9)", "18:2(n-6)", 
    "20:4(n-6)", "22:5(n-6)", "18:3(n-3)", "18:4(n-3)", "20:5(n-3)", 
    "22:5(n-3)", "22:6(n-3)")

perch_FA_prop_PCA_species$FA <- trimmed.FA.names

perch_FA_prop_PCA_site <- cbind(perch_FA_prop_covariates,
                                perch_FA_prop_PCA_site[, 1:2])

#### Plot ----

png("Perch FA Proportions PCA.png",
    width = 19,
    height= 12, 
    units = "cm",
    res = 500)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_point(data = perch_FA_prop_PCA_site,
             aes(x = PC1,
                 y = PC2,
                 fill = as.factor(MPconcentration)),
             size = 2,
             alpha = 0.75,
             shape = 21) +
  geom_text(data = perch_FA_prop_PCA_species,
            aes(x = PC1, 
                y = PC2, 
                label = FA),
            alpha = 0.95,
            size = 7 / .pt,
            colour = "black") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "inferno") +
  labs(x = "PC1", 
       y = "PC2") +
  theme1 +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()

### RDA ----

perch_FA_prop_rda <- 
  rda(trimmed_perch_FA_clr ~ log(MPconcentration + 1) + body.weight,
      scale. = FALSE,
      data = perch_FA_prop_covariates)

summary(perch_FA_prop_rda)

anova(perch_FA_prop_rda, by = "term")
anova(perch_FA_prop_rda, by = "margin")
anova(perch_FA_prop_rda, by = "onedf")

### CCA ----

# Add diet data

perch_FA_prop_covariates <- 
  left_join(perch_FA_prop_covariates, 
            diet_summary, 
            by = "corral")

perch_FA_prop_covariates$corral <- as.factor(perch_FA_prop_covariates$corral)

perch_FA_prop_cca <- 
  cca((trimmed_perch_FA) ~ 
        corral + 
        body.weight,
      data = perch_FA_prop_covariates)

summary(perch_FA_prop_cca)

RsquareAdj(perch_FA_prop_cca)

screeplot(perch_FA_prop_cca)

anova(perch_FA_prop_cca, by = "term")
anova(perch_FA_prop_cca, by = "margin")
anova(perch_FA_prop_cca, by = "onedf")

plot(perch_FA_prop_cca, 
     scaling = 1, 
     display = c("lc", "cn"))

plot(perch_FA_prop_cca, 
     scaling = 2, 
     display = c("sp","cn"))

summary(perch_FA_prop_cca)

#### Pull out scores for scaling 1----

# 'Site' scores
perch_FA_prop_cca_site1 <- 
  as.data.frame(scores(perch_FA_prop_cca, display = "site",
                       scaling = 1))

# 'Species' scores
perch_FA_prop_cca_species1 <- 
  as.data.frame(scores(perch_FA_prop_cca, display = "species",
                       scaling = 1))

trimmed.FA.names <- 
  c("14:0", "16:0", "18:0", "16:1(n-7)", "18:(1n-7)", "18:1(n-9)", "18:2(n-6)", 
    "20:4(n-6)", "22:5(n-6)", "18:3(n-3)", "18:4(n-3)", "20:5(n-3)", 
    "22:5(n-3)", "22:6(n-3)")

perch_FA_prop_cca_species1$FA <- trimmed.FA.names

perch_FA_prop_cca_site1 <- cbind(perch_FA_prop_covariates,
                                 perch_FA_prop_cca_site1[, 1:2])

perch_FA_prop_cca_centroids1 <- 
  as.data.frame(scores(perch_FA_prop_cca, 
                       display = "cn",
                       scaling = 1))

perch_FA_prop_cca_centroids1$vars <- c("0(B)",
                                      "414(C)",
                                      "29,240(D)",
                                      "100(E)",
                                      "6(F)",
                                      "7,071(G)",
                                      "0(H)",
                                      "1,710(I)")

perch_FA_prop_cca_centroids1$vars <-
  factor(perch_FA_prop_cca_centroids1$vars,
         levels = c("0(B)",
                    "0(H)",
                    "6(F)",
                    "100(E)",
                    "414(C)",
                    "1,710(I)",
                    "7,071(G)",
                    "29,240(D)"))

perch_FA_prop_cca_centroids1$corral <- as.factor(LETTERS[2:9])
perch_FA_prop_cca_centroids1 <-
  rename(perch_FA_prop_cca_centroids1,
         cCCA1 = CCA1,
         cCCA2 = CCA2)

perch_FA_prop_cca_site1 <-
  left_join(perch_FA_prop_cca_site1,
            perch_FA_prop_cca_centroids1,
            by = "corral")

perch_FA_prop_cca_bp1 <- 
  as.data.frame(scores(perch_FA_prop_cca, 
                       display = "bp",
                       scaling = 1))

perch_FA_prop_cca_bp1 <- perch_FA_prop_cca_bp1[8,]

perch_FA_prop_cca_bp1$label <- "Body Weight"

#### Pull out scores for scaling 2----

# 'Site' scores
perch_FA_prop_cca_site2 <- 
  as.data.frame(scores(perch_FA_prop_cca, display = "site",
                       scaling = 2))

# 'Species' scores
perch_FA_prop_cca_species2 <- 
  as.data.frame(scores(perch_FA_prop_cca, display = "species",
                       scaling = 2))

perch_FA_prop_cca_species2$FA <- trimmed.FA.names

perch_FA_prop_cca_site2 <- cbind(perch_FA_prop_covariates,
                                 perch_FA_prop_cca_site2[, 1:2])

perch_FA_prop_cca_centroids2 <- 
  as.data.frame(scores(perch_FA_prop_cca, display = "cn",
                       scaling = 2))

perch_FA_prop_cca_centroids2$vars <- c("0(B)",
                                       "414(C)",
                                       "29,240(D)",
                                       "100(E)",
                                       "6(F)",
                                       "7,071(G)",
                                       "0(H)",
                                       "1,710(I)")

perch_FA_prop_cca_centroids2$vars <-
  factor(perch_FA_prop_cca_centroids2$vars,
         levels = c("0(B)",
                    "0(H)",
                    "6(F)",
                    "100(E)",
                    "414(C)",
                    "1,710(I)",
                    "7,071(G)",
                    "29,240(D)"))

perch_FA_prop_cca_centroids2$corral <- as.factor(LETTERS[2:9])
perch_FA_prop_cca_centroids2 <-
  rename(perch_FA_prop_cca_centroids2,
         cCCA1 = CCA1,
         cCCA2 = CCA2)

perch_FA_prop_cca_site2 <-
  left_join(perch_FA_prop_cca_site2,
            perch_FA_prop_cca_centroids2,
            by = "corral")

perch_FA_prop_cca_bp2 <- 
  as.data.frame(scores(perch_FA_prop_cca, 
                       display = "bp",
                       scaling = 2))

perch_FA_prop_cca_bp2 <- perch_FA_prop_cca_bp2[8,]

perch_FA_prop_cca_bp2$label <- "Body Weight"

#### Plot scaling 1 ----

png("Perch FA Proportions CCA Spider Scaling 1.png",
    width = 18,
    height= 10, 
    units = "cm",
    res = 600)
  
ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_segment(data = perch_FA_prop_cca_site1,
             aes(x = cCCA1,
                 y = cCCA2,
                 xend = CCA1,
                 yend = CCA2,
                 colour = vars),
             alpha = 0.75,
             linewidth = 0.5,
             arrow = arrow(angle = 20,
                           length = unit(0.25, "cm"),
                           type = "open")) +
  geom_point(data = perch_FA_prop_cca_site1,
               aes(x = cCCA1,
                   y = cCCA2,
                   fill = vars),
               size = 2,
             shape = 21,
             alpha = 0.95,
             colour = "grey30") +
  geom_text(data = perch_FA_prop_cca_bp1,
             aes(x = CCA1*4,
                 y = CCA2*20,
                 label = label),
             size = 10 / .pt,
             colour =  "blue3") +
  geom_segment(data = perch_FA_prop_cca_bp1,
               aes(x = 0,
                   y = 0,
                   xend = CCA1*4,
                   yend = CCA2*4,
                   colour = vars),
               alpha = 0.75,
               linewidth = 0.5,
               arrow = arrow(angle = 20,
                             length = unit(0.25, "cm"),
                             type = "open"),
               colour = "blue3") +
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
  theme(legend.key.size = unit(0.2, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()


#### Plot scaling 2 ----

png("Perch FA Proportions CCA Spider Scaling 2.png",
    width = 18,
    height= 10, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_point(data = perch_FA_prop_cca_site2,
             aes(x = cCCA1,
                 y = cCCA2,
                 fill = vars),
             size = 2,
             shape = 21,
             alpha = 0.95,
             colour = "grey30") +
  geom_text(data = perch_FA_prop_cca_species2,
            aes(x = CCA1,
                y = CCA2,
                label = FA),
            size = 10 / .pt,
            colour = "blue3") +
  geom_text(data = perch_FA_prop_cca_bp2,
            aes(x = CCA1,
                y = CCA2*10,
                label = label),
            size = 10 / .pt,
            colour =  "blue3") +
  geom_segment(data = perch_FA_prop_cca_bp2,
               aes(x = 0,
                   y = 0,
                   xend = CCA1,
                   yend = CCA2,
                   colour = vars),
               alpha = 0.75,
               linewidth = 0.5,
               arrow = arrow(angle = 20,
                             length = unit(0.25, "cm"),
                             type = "open"),
               colour = "blue3") +
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
  theme(legend.key.size = unit(0.2, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()


### nMDS ----

# Calculate distance matrix
perch_FA_prop_diss <- as.matrix(vegdist(trimmed_perch_FA, 
                                   method = "bray", 
                                   na.rm = TRUE), 
                           labels = TRUE)

# NMDS_scree(perch_FA_prop_diss)  # 3 dimensions looks good

set.seed(5465)

perch_FA_prop_nMDS1 <- 
  metaMDS(perch_FA_prop_diss,
          distance = "bray",
          k = 3,
          trymax = 250,
          wascores = TRUE,
          expand = TRUE,
          autotransform = FALSE)

# Shepards test/goodness of fit
goodness(perch_FA_prop_nMDS1)
stressplot(perch_FA_prop_nMDS1)

perch_FA_prop_scores <- `sppscores<-`(perch_FA_prop_nMDS1, 
                                      trimmed_perch_FA)

perch_FA_prop_data.scores <- as.data.frame(perch_FA_prop_scores$points)

perch_FA_prop_data.scores2 <- cbind(perch_FA_prop_data.scores,
                                    perch_FA_prop_covariates)

perch_FA_prop_variable_scores <- 
  as.data.frame(perch_FA_prop_scores$species)

perch_FA_prop_variable_scores$FA <- 
  trimmed.FA.names

#### Generate hulls ----

perch_FA_prop_data.scores2$MPconcentration <- 
  as.factor(perch_FA_prop_data.scores2$MPconcentration)

perch_FA_prop_hulls <- data.frame()

for(i in 1:length(unique(perch_FA_prop_data.scores2$MPconcentration))) {
  subset <- 
    perch_FA_prop_data.scores2[perch_FA_prop_data.scores2$MPconcentration ==
                                 unique(perch_FA_prop_data.scores2$MPconcentration)[i],
                               ]
  hullrows <- chull(subset[,c(1:2)])
  hull <- subset[hullrows,]
  perch_FA_prop_hulls <- rbind(perch_FA_prop_hulls, hull)
}

#### Plot ----

pal2 <-
  colorRampPalette(c("#EBEA9A", "#E4C008", "#A17244", "#42511A", 
                     "#1B334D", "#151918"))(7)

png("Perch FA Proportions nMDS Plot.png",
    width = 12,
    height= 12, 
    units = "cm",
    res = 500)

ggplot() +
  geom_polygon(data = perch_FA_prop_hulls,
               aes(x = MDS1,
                   y = MDS2,
                   fill = MPconcentration,
                   colour = MPconcentration),
               alpha = 0.5,
               size = 0.5) +
  geom_segment(data = perch_FA_prop_variable_scores,
               aes(x = 0, y = 0, xend = MDS1*0.9, yend = MDS2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple3",
               alpha = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = perch_FA_prop_hulls,
             aes(x = MDS1,
                 y = MDS2,
                 fill = MPconcentration),
             size = 2,
             alpha = 0.75,
             shape = 21) +
  geom_text(data = perch_FA_prop_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = FA),
            alpha = 0.9,
            size = 7 / .pt,
            colour = "purple3") +
  scale_fill_manual(values = pal2,
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_manual(values = pal2,
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1 +
  theme(legend.position = "bottom")

dev.off()

### PERMANOVA ####

perm <- how(nperm = 999)
setBlocks(perm) <- perch_FA_prop_covariates$corral

perch.permanova <- 
  adonis2(perch_FA_prop_diss ~ log(MPconcentration + 1) + body.weight,
          permutations = perm,
          data = perch_FA_prop_covariates,
          method = "bray",
          parallel = 8)

perch.permanova  

perch.betadisper <- 
  betadisper(vegdist(trimmed_perch_FA, 
                     method = "bray", 
                     na.rm = TRUE),
             group = as.factor(perch_FA_prop_covariates$MPconcentration),
             type = "median")

anova(perch.betadisper)
boxplot(perch.betadisper)
TukeyHSD(perch.betadisper)

# No significant effect of MP concentration or body weight

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
  filter(FA == "C_18.1n.9" |
           FA == "C_18.2n.6" |
           FA == "C_18.3n.3" |
           FA == "C_18.4n.3" |
           FA == "C_20.4n.6" |
           FA == "C_22.5n.6" |
           FA == "C_22.6n.3" )

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
  facet_wrap(~ FA, scale = "free") +
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
            as.factor(date), 
          data = zoop.end)

plot(simulateResiduals(DHA.prop.zoop.mod1))

summary(DHA.prop.zoop.mod1)  
# no effect of MP concentration
# more DHA at the end

#### EPA ----

EPA.prop.zoop.mod1 <- 
  glmmTMB(C_20.5n.3 ~ 
            log(MPconcentration + 6) +
            as.factor(date), 
          data = zoop.end)

plot(simulateResiduals(EPA.prop.zoop.mod1))

summary(EPA.prop.zoop.mod1)  
# no effect of MP concentration
# less EPA at the end

#### ARA ----

ARA.prop.zoop.mod1 <- 
  glmmTMB(C_20.4n.6 ~ 
            MPconcentration +
            as.factor(date), 
          data = zoop.end)

plot(simulateResiduals(ARA.prop.zoop.mod1))

summary(ARA.prop.zoop.mod1)  
# no effect of MP concentration
# less ARA at the end

zoop.end.long <-
  zoop.end %>% 
  select(-c(7:30,32:39,41,43:50)) %>% 
  pivot_longer(cols = c(7:9),
               names_to = "FA")
png("Zoop Compositional HUFA Plot.png",
    width = 8.84,
    height= 6, 
    units = "cm",
    res = 500)

ggplot(zoop.end.long) +
  geom_boxplot(aes(x = as.factor(date),
                  y = value,
                  fill = FA)) +
  scale_fill_viridis_d(name = "",
                       labels = c("ARA",
                                  "EPA",
                                  "DHA")) +
  labs(x = "Date",
       y = "Proportion") +
  theme1

dev.off()

### PCA ----

#### Pull out covariates ----

# Remove any FAs that contribute less than 1% and re-scale response so it 
# sums to 1

trimmed_zoop_FA <- 
  data.frame(zoop_FA_prop[,c(7:16,18:25,27:34,36:44)]) %>% 
  select(where(function(x){mean(x) >= 0.01}))

# trimmed_zoop_FA <- trimmed_zoop_FA/rowSums(trimmed_zoop_FA)

zoop_FA_prop_covariates <- zoop_FA_prop[,c(1:3,51:54),]

zoop_FA_prop_totals <- zoop_FA_prop[,c(17,26,35,45)]

summary(zoop_FA_prop_totals)

# Convert to CLRs

trimmed_zoop_FA_clr <- decostand(trimmed_zoop_FA,
                                 method = "clr",
                                 MARGIN = 1)



#### Run PCA ----
zoop_FA_prop_pca <- rda(trimmed_zoop_FA_clr,
                        scale. = FALSE)

# Bar plot of relative eigenvalues
barplot(as.vector(zoop_FA_prop_pca$CA$eig)/sum(zoop_FA_prop_pca$CA$eig))

# Calculate percentage of variance explained by first 2 aaxes
sum((as.vector(zoop_FA_prop_pca$CA$eig)/sum(zoop_FA_prop_pca$CA$eig))[1:2])

summary(zoop_FA_prop_pca)

# 'Site' scores
zoop_FA_prop_PCA_site <- zoop_FA_prop_pca$CA$u

# 'Species' scores
zoop_FA_prop_PCA_species <- data.frame(zoop_FA_prop_pca$CA$v)

trimmed.FA.names.zoops <- 
  c("14:0", "16:0", "18:0", "16:1(n-7)", "18:(1n-7)", "18:1(n-9)", "22:1(n-9)",
    "18:2(n-6)", "18:3(n-6)", "20:4(n-6)", "22:5(n-6)", "18:3(n-3)", 
    "18:4(n-3)", "20:5(n-3)", "22:6(n-3)")

zoop_FA_prop_PCA_species$FA <- trimmed.FA.names.zoops

zoop_FA_prop_PCA_site <- cbind(zoop_FA_prop_covariates,
                               zoop_FA_prop_PCA_site[, 1:2])

#### Plot ----

png("Zooplankton FA Proportions PCA.png",
    width = 19,
    height= 12, 
    units = "cm",
    res = 500)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_point(data = zoop_FA_prop_PCA_site,
             aes(x = PC1,
                 y = PC2,
                 fill = as.factor(MPconcentration),
                 shape = as.factor(date)),
             size = 2,
             alpha = 0.75) +
  geom_text(data = zoop_FA_prop_PCA_species,
            aes(x = PC1, 
                y = PC2, 
                label = FA),
            alpha = 0.95,
            size = 7 / .pt,
            colour = "black") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "inferno") +
  scale_shape_manual(name = "Date",
                     values = c(21:23)) +
  labs(x = "PC1",
       y = "PC2") +
  theme1 +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()

### RDA ----

zoop_FA_prop_rda <- 
  rda(trimmed_zoop_FA_clr ~ as.factor(MPconcentration) + as.factor(date),
      scale. = FALSE,
      data = zoop_FA_prop_covariates)

summary(zoop_FA_prop_rda)

zoop_FA_prop_rda2 <- 
  rda(trimmed_zoop_FA_clr ~ as.factor(MPconcentration),
      scale. = FALSE,
      data = zoop_FA_prop_covariates)

zoop_FA_prop_rda3 <- 
  rda(trimmed_zoop_FA_clr ~ as.factor(date),
      scale. = FALSE,
      data = zoop_FA_prop_covariates)

anova(zoop_FA_prop_rda, by = "term")
anova(zoop_FA_prop_rda, by = "margin")
anova(zoop_FA_prop_rda, by = "onedf")

### CCA ----

zoop_FA_prop_cca <- 
  cca(trimmed_zoop_FA ~ corral + as.factor(date),
      scale. = FALSE,
      data = zoop_FA_prop_covariates)

plot(zoop_FA_prop_cca, scaling = 1, display = c("lc", "cn"))

plot(zoop_FA_prop_cca, scaling = 2, display = c("sp", "cn"))

screeplot(zoop_FA_prop_cca)

summary(zoop_FA_prop_cca)

set.seed(4399)

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

trimmed.FA.names.zoops <- 
  c("14:0", "16:0", "18:0", 
    "16:1(n-7)", "18:1(n-7)", "18:1(n-9)", 
    "22:1(n-9)",  "18:2(n-6)", "18:3(n-6)", 
    "20:4(n-6)", "22:5(n-6)", "18:3(n-3)", 
    "18:4(n-3)", "20:5(n-3)", "22:6(n-3)")

zoop_FA_prop_cca_species1$FA <- trimmed.FA.names.zoops

zoop_FA_prop_cca_site1 <- cbind(zoop_FA_prop_covariates,
                                zoop_FA_prop_cca_site1[, 1:2])

zoop_FA_prop_cca_site1$label <- zoop_FA_prop_cca_site1$corral

levels(zoop_FA_prop_cca_site1$label) <- c("24(A)",
                                         "0(B)",
                                         "414(C)",
                                         "29,240(D)",
                                         "100(E)",
                                         "6(F)",
                                         "7,071(G)",
                                         "0(H)",
                                         "1,710(I)")

zoop_FA_prop_cca_site1$label <-
  factor(zoop_FA_prop_cca_site1$label,
         levels = c("0(B)",
                    "0(H)",
                    "6(F)",
                    "24(A)",
                    "100(E)",
                    "414(C)",
                    "1,710(I)",
                    "7,071(G)",
                    "29,240(D)"))


zoop_FA_prop_cca_site1$point <- 
  as.factor(as.character(zoop_FA_prop_cca_site1$date))

levels(zoop_FA_prop_cca_site1$point) <- c("Start",
                                         "Mid-point",
                                         "End")
# Generate hulls
zoop_FA_date_hulls <- data.frame()

for(i in 1:length(unique(zoop_FA_prop_cca_site1$point))) {
  subset <- 
    zoop_FA_prop_cca_site1[zoop_FA_prop_cca_site1$point ==
                            unique(zoop_FA_prop_cca_site1$point)[i],
    ]
  hullrows <- chull(subset[,c(8:9)])
  hull <- subset[hullrows,]
  zoop_FA_date_hulls <- rbind(zoop_FA_date_hulls, hull)
}

zoop_FA_prop_cca_cn1 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 1))

zoop_FA_prop_cca_cn1$label <-
  factor(c("24(A)",
              "0(B)",
              "414(C)",
              "29, 240(D)",
              "100(E)",
              "6(F)",
              "7, 071(G)",
              "0(H)",
              "1, 710(I)",
              "Start",
              "Mid-point",
    "End"))

##### Plot ----

png("Zooplankton FA Scaling 1 Proportions CCA.png",
    width = 19,
    height= 8, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_polygon(data = zoop_FA_date_hulls,
               aes(x = CCA1,
                   y = CCA2,
                   group = point),
               alpha = 0.25,
               linetype = "dashed",
               linewidth = 0.5) +
  geom_point(data = zoop_FA_prop_cca_site1,
             aes(x = CCA1,
                 y = CCA2,
                 fill = label,
                 shape = point),
             size = 2,
             alpha = 0.75,
             colour = "grey30") +
  geom_text_repel(data = zoop_FA_prop_cca_cn1,
                  aes(x = CCA1,
                      y = CCA2,
                      label = label),
                  size = 10 / .pt,
                  colour =  "blue3",
                  box.padding = unit(0, "cm")) +
  scale_shape_manual(values = c(21,24,22),
                     name = "Experimental Time Point") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  labs(x = "CCA1", 
       y = "CCA2") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme1 +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()

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

levels(zoop_FA_prop_cca_site1$label) <- c("24(A)",
                                          "0(B)",
                                          "414(C)",
                                          "29,240(D)",
                                          "100(E)",
                                          "6(F)",
                                          "7,071(G)",
                                          "0(H)",
                                          "1,710(I)")

zoop_FA_prop_cca_site2$label <-
  factor(zoop_FA_prop_cca_site2$label,
         levels = c("0(B)",
                    "0(H)",
                    "6(F)",
                    "24(A)",
                    "100(E)",
                    "414(C)",
                    "1,710(I)",
                    "7,071(G)",
                    "29,240(D)"))


zoop_FA_prop_cca_site2$point <- 
  as.factor(as.character(zoop_FA_prop_cca_site2$date))

levels(zoop_FA_prop_cca_site2$point) <- c("Start",
                                          "Mid-point",
                                          "End")

zoop_FA_prop_cca_cn1 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 1))

zoop_FA_prop_cca_cn1$label <-
  factor(c("24(A)",
              "0(B)",
              "414(C)",
              "29, 240(D)",
              "100(E)",
              "6(F)",
              "7, 071(G)",
              "0(H)",
              "1, 710(I)",
              "Start",
              "Mid-point",
    "End"))

zoop_FA_prop_cca_cn2 <-
  data.frame(scores(zoop_FA_prop_cca,
                    display = "cn",
                    scaling = 2))

zoop_FA_prop_cca_cn2$label <-
  factor(c("24(A)",
           "0(B)",
           "414(C)",
           "29, 240(D)",
           "100(E)",
           "6(F)",
           "7, 071(G)",
           "0(H)",
           "1, 710(I)",
           "Start",
           "Mid-point",
           "End"))



##### Plot ----

png("Zooplankton FA Scaling 2 Proportions CCA.png",
    width = 12,
    height= 8, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             linewidth = 0.25) +
  geom_text_repel(data = zoop_FA_prop_cca_species2,
                  aes(x = CCA1,
                      y = CCA2,
                      label = FA),
                  size = 8 / .pt,
                  colour = "blue3",
                  box.padding = 0) +
  geom_text_repel(data = zoop_FA_prop_cca_cn2,
                  aes(x = CCA1,
                      y = CCA2,
                      label = label),
                  size = 10 / .pt,
                  colour =  "purple3",
                  box.padding = unit(0, "cm")) +
  scale_shape_manual(values = c(21,24,22),
                     name = "Experimental Time Point") +
  scale_fill_viridis_d(name =
                         expression(paste("Exposure Concentration (MPs" ~
                                            L ^ -1 * ")")),
                       option = "plasma") +
  scale_y_continuous(limits = c(-1.25, 1)) +
  labs(x = "CCA1", 
       y = "CCA2") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme1 +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0, "cm"))

dev.off()



### nMDS ----

# Calculate distance matrix
zoop_FA_prop_diss <- as.matrix(vegdist(trimmed_zoop_FA, 
                                  method = "bray", 
                                  na.rm = TRUE), 
                          labels = TRUE)

set.seed(5465)

zoop_FA_prop_nMDS1 <- 
  metaMDS(zoop_FA_prop_diss,
          distance = "bray",
          k = 4,
          trymax = 250,
          wascores = TRUE,
          expand = TRUE,
          autotransform = FALSE)

# Shepards test/goodness of fit
goodness(zoop_FA_prop_nMDS1)
stressplot(zoop_FA_prop_nMDS1)

zoop_FA_prop_scores <- `sppscores<-`(zoop_FA_prop_nMDS1, trimmed_zoop_FA)

zoop_FA_prop_data.scores <- as.data.frame(zoop_FA_prop_scores$points)

zoop_FA_prop_data.scores2 <- cbind(zoop_FA_prop_data.scores,
                                   zoop_FA_prop_covariates)

zoop_FA_variable_prop_scores <- 
  as.data.frame(zoop_FA_prop_scores$species)

zoop_FA_variable_prop_scores$FA <- 
  trimmed.FA.names.zoops

#### Generate hulls ----

zoop_FA_prop_data.scores2$MPconcentration <- 
  as.factor(zoop_FA_prop_data.scores2$MPconcentration)

zoop_FA_prop_data.scores2$date <- 
  as.factor(zoop_FA_prop_data.scores2$date)

zoop_FA_prop_hulls <- data.frame()

for(i in 1:length(unique(zoop_FA_prop_data.scores2$date))) {
  hull <-
    zoop_FA_prop_data.scores2[zoop_FA_prop_data.scores2$date ==
                           unique(zoop_FA_prop_data.scores2$date)[i],
    ][chull(zoop_FA_prop_data.scores2[zoop_FA_prop_data.scores2$date ==
                                   unique(zoop_FA_prop_data.scores2$date)[i],
                                 c(1:2)]),]
  zoop_FA_prop_hulls <- rbind(zoop_FA_prop_hulls, hull)
}

#### Plot ----

png("Zooplankton Proportions MDS Plot.png",
    width = 19,
    height= 12, 
    units = "cm",
    res = 500)

ggplot() +
  geom_polygon(data = zoop_FA_prop_hulls,
               aes(x = MDS1,
                   y = MDS2,
                   colour = date),
               alpha = 0.1,
               size = 1) +
  geom_segment(data = zoop_FA_variable_prop_scores,
               aes(x = 0, y = 0, xend = MDS1*0.9, yend = MDS2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple4",
               alpha = 0.75) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = zoop_FA_prop_hulls,
             aes(x = MDS1,
                 y = MDS2,
                 shape = date,
                 fill = MPconcentration),
             alpha = 0.75,
             size = 2) +
  geom_text(data = zoop_FA_variable_prop_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = FA),
            alpha = 0.9,
            size = 7 / .pt,
            colour = "purple4") +
  scale_fill_brewer(type = "seq",
                    palette = "YlOrRd",
                    name =
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "qual",
                    palette = "Set3",
                      name = "Date") +
  scale_shape_manual(name = "Date",
                     values = c(21:23)) +
  theme1

dev.off()

### PERMANOVA ####

perm <- how(nperm = 999)
setBlocks(perm) <- zoop_FA_prop_covariates$corral

zoop.permanova <- 
  adonis2(zoop_FA_prop_diss ~ log(MPconcentration + 1) * scaled.date,
          permutations = perm,
          data = zoop_FA_prop_covariates,
          method = "bray",
          parallel = 8)

zoop.permanova

zoop.betadisper <- 
  betadisper(vegdist(trimmed_zoop_FA, 
                     method = "bray", 
                     na.rm = TRUE),
             group = as.factor(zoop_FA_prop_covariates$MPconcentration),
             type = "centroid")

anova(zoop.betadisper)
plot(zoop.betadisper)
boxplot(zoop.betadisper)
TukeyHSD(zoop.betadisper)

# Plot proportions of major groups ----
## Fish food ----

food_FA_prop_long <- 
  food_FA_prop %>% 
  pivot_longer(cols = c(17,26,35,45),
               names_to = "FA")

food_FA_prop_long$FA <- as.factor(food_FA_prop_long$FA)

levels(food_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

png("Fish Food FAs Proportionss Plot.png",
    width = 8.84,
    height= 4, 
    units = "cm",
    res = 600)

ggplot(food_FA_prop_long) +
  geom_col(aes(x = ID,
               y = value,
               fill = FA)) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  labs(x = "Sample",
       y = "Proportion") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

dev.off()

## Perch ----

perch_FA_prop_long <- 
  perch_FA_prop2 %>% 
  pivot_longer(cols = c(17,26,35,45),
               names_to = "FA")

perch_FA_prop_long$FA <- as.factor(perch_FA_prop_long$FA)

levels(perch_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

png("Perch FAs Proportionss Plot.png",
    width = 18,
    height= 8.84, 
    units = "cm",
    res = 600)

ggplot(perch_FA_prop_long) +
  geom_col(aes(x = ID,
               y = value,
               fill = FA)) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  labs(x = "Sample",
       y = "Proportion") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  facet_wrap(~MPconcentration,
             scales = "free_x",
             nrow = 2) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

dev.off()

## Zooplankton ----

zoop_FA_prop_long <- 
  zoop_FA_prop %>% 
  pivot_longer(cols = c(17,26,35,45),
               names_to = "FA")

zoop_FA_prop_long$FA <- as.factor(zoop_FA_prop_long$FA)

levels(zoop_FA_prop_long$FA) <-
  c("MUFA",
    "n-3 PUFA",
    "n-6 PUFA",
    "SFA")

png("Zooplankton FAs Proportionss Plot.png",
    width = 18,
    height= 12, 
    units = "cm",
    res = 600)

ggplot(zoop_FA_prop_long) +
  geom_col(aes(x = reorder(ID, date),
               y = value,
               fill = FA,
               colour = as.factor(date)),
           size = 0.75) +
  scale_fill_viridis_d(option = "turbo",
                       name = "") +
  scale_colour_viridis_d(name = "Date",
                         option = "plasma") +
  labs(x = "Sample",
       y = "Proportion") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  facet_wrap(~ MPconcentration,
             scales = "free_x") +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

dev.off()

