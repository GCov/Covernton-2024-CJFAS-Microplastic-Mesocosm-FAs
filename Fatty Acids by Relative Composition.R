# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(plyr)
library(dplyr)
library(MuMIn)
library(vegan)
library(tidyr)
library(Hmsc)

# Load data ----

FAs_percent <- read.csv("FAs_percent.csv", header = TRUE)
str(FAs_percent)
FAs_percent$ID <- as.factor(FAs_percent$ID)
FAs_percent$corral <- as.factor(FAs_percent$corral)
FAs_percent$sample.type <- as.factor(FAs_percent$sample.type)
FAs_percent$date <- as.Date(FAs_percent$date,
                    format = "%b. %d, %Y")

FAs_percent <- left_join(FAs_percent, 
                 treatments,
                 by = "corral")

# Assign corral A MP concentration of 24 for zoop_FA

FAs_percent$MPconcentration[is.na(FAs_percent$MPconcentration)] <- 24

# Convert to proportions
FAs_prop <- FAs_percent
FAs_prop[7:50] <- FAs_percent[7:50] / 100

# Separate out by sample type ----

perch_FA_prop <- subset(FAs_prop, sample.type == "perch")
zoop_FA_prop <- subset(FAs_prop, sample.type == "zooplankton")
food_FA_prop <- subset(FAs_prop, sample.type == "fish food")

# Add in perch biometrics data ----

# Combine with perch FA data

perch_FA_prop2 <- left_join(perch_FA_prop, 
                               perch_biometrics, 
                       by = c("ID", "corral"))

# Add in perch population data ----

perch_pop <- read.csv("fish_pop.csv", header = TRUE)

# Combine with perch FA data

perch_FA_prop2 <- left_join(perch_FA_prop2,
                            perch_pop,
                            by = c("corral", "MPconcentration"))

# Scale and center date

perch_FA_prop2$date2 <- 
  as.numeric(scale(as.numeric(perch_FA_prop2$date), center = TRUE))

## Plot by fatty acid composition ----

# Put data into long form and remove totals

perch_FA_prop_long <- 
  perch_FA_prop2[,c(1:3,6,7:16,18:25,27:34,36:44,51:59)] %>%
  pivot_longer(names(perch_FA_prop2[c(7:16,18:25,27:34,36:44)]),
               names_to = "fatty.acid",
               values_to = "prop")

# Plot

ggplot(perch_FA_prop_long) +
  geom_col(aes(x = ID,
               y = prop,
               fill = fatty.acid),
           colour = "black",
           size = 0.25) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1)) +
  labs(x = expression(paste("MP Exposure Concentration (particles "~L^-1*")")),
       y = "prop Composition") +
  scale_fill_hue(name = "Fatty Acid",
                 l = 75) +
  facet_grid(.~MPconcentration,
             scales = "free_x",
             switch = "x",
             space = "free_x") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Explore different indicators ----

# Scale and center zooplankton date

zoop_FA_prop$date2 <- 
  as.numeric(scale(as.numeric(zoop_FA_prop$date), center = TRUE))

## Zooplankton PUFAs ----

zoop_FA_prop_PUFA_mod1 <- 
  glmmTMB(PUFAs ~ log(MPconcentration + 1),
          data = subset(zoop_FA_prop, date > "2021-08-01"),
          family = beta_family(link = "logit"))

plotResiduals(zoop_FA_prop_PUFA_mod1)

summary(zoop_FA_prop_PUFA_mod1)

## Zooplankton n-6 PUFAs ----

zoop_FA_prop_n.6.PUFA_mod1 <- 
  glmmTMB(total_N.6_PUFAs ~ log(MPconcentration + 1),
          data = subset(zoop_FA_prop, date > "2021-08-01"),
          family = beta_family(link = "logit"))

plotResiduals(zoop_FA_prop_n.6.PUFA_mod1)

summary(zoop_FA_prop_n.6.PUFA_mod1)

## Zooplankton n-3 PUFAs ----

zoop_FA_prop_n.3.PUFA_mod1 <- 
  glmmTMB(total_N.3_PUFAs ~ log(MPconcentration + 1),
          data = subset(zoop_FA_prop, date > "2021-08-01"),
          family = beta_family(link = "logit"))

plotResiduals(zoop_FA_prop_n.3.PUFA_mod1)

summary(zoop_FA_prop_n.3.PUFA_mod1)

# Zooplankton nMDS ----

# Pull out covariates

zoop_FA_covariates <- zoop_FA_prop[,c(1:3,51),]

zoop_FA_totals <- zoop_FA_prop[,c(17,26,35,45)]

# Pull out fatty acid composition matrix

zoop_FA_prop_matrix <- zoop_FA_prop[,c(7:16,18:25,27:34,36:44)]

# Calculate distance matrix
zoop_FA_diss <- as.matrix(vegdist(zoop_FA_prop_matrix, 
                                   method = "euclidean", 
                                   na.rm = TRUE), 
                           labels = TRUE)

NMDS_scree <-
  function(x) {
    #where x is the name of the data frame variable
    plot(
      rep(1, 10),
      replicate(10, metaMDS(
        x, autotransform = F, k = 1
      )$stress),
      xlim = c(1, 10),
      ylim = c(0, 0.30),
      xlab = "# of Dimensions",
      ylab = "Stress",
      main = "NMDS stress plot"
    )
    for (i in 1:10) {
      points(rep(i + 1, 10), replicate(10, metaMDS(
        x, autotransform = F, k = i + 1
      )$stress))
    }
  }

NMDS_scree(zoop_FA_diss)  # 4 dimensions looks good

set.seed(5465)

zoop_FA_nMDS1 <- 
  metaMDS(zoop_FA_diss,
          distance = "euclidean",
          k = 4,
          trymax = 250,
          wascores = TRUE,
          expand = TRUE,
          autotransform = FALSE)

# Shepards test/goodness of fit
goodness(zoop_FA_nMDS1)
stressplot(zoop_FA_nMDS1)

zoop_FA_data.scores <- as.data.frame(scores(zoop_FA_nMDS1))
zoop_FA_data.scores2 <- cbind(zoop_FA_data.scores,
                               zoop_FA_covariates)

zoop_FA_scores <- `sppscores<-`(zoop_FA_nMDS1, zoop_FA_prop_matrix)

zoop_FA_variable_scores <- 
  as.data.frame(zoop_FA_scores$species)

zoop_FA_variable_scores$variable <- 
  rownames(zoop_FA_variable_scores)

## Generate hulls ----

zoop_FA_data.scores2$MPconcentration <- 
  as.factor(zoop_FA_data.scores2$MPconcentration)

zoop_FA_hulls <- data.frame()

for(i in 1:length(unique(zoop_FA_data.scores2$MPconcentration))) {
  hull <-
    zoop_FA_data.scores2[zoop_FA_data.scores2$MPconcentration ==
                            unique(zoop_FA_data.scores2$MPconcentration)[i],
    ][chull(zoop_FA_data.scores2[zoop_FA_data.scores2$MPconcentration ==
                                    unique(zoop_FA_data.scores2$MPconcentration)[i],
                                  c(1:2)]),]
  zoop_FA_hulls <- rbind(zoop_FA_hulls, hull)
}

## Plot ----

png("Zooplankton MDS Plot.png",
    width = 19,
    height= 15, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = zoop_FA_hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   fill = MPconcentration,
                   colour = MPconcentration),
               alpha = 0.75,
               size = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_label(data = zoop_FA_hulls,
             aes(x = NMDS1,
                 y = NMDS2,
                 label = ID),
             alpha = 0.3,
             size = 4) +
  geom_text(data = zoop_FA_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = variable,
                angle = MDS1*MDS2*0.6),
            alpha = 0.9,
            size = 3,
            colour = "purple") +
  scale_fill_brewer(type = "seq",
                    palette = "YlGnBu",
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "seq",
                      palette = "YlGnBu",
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1

dev.off()

# PCA

zoop_FA_pca <- prcomp(zoop_FA_prop_matrix, 
                       center = TRUE,
                       scale. = TRUE)

summary(zoop_FA_pca)

library(ggbiplot)

ggbiplot(zoop_FA_pca,
         labels = zoop_FA_covariates$ID,
         ellipse = TRUE,
         groups = as.factor(zoop_FA_covariates$MPconcentration))


# Perch nMDS ----

# Pull out covariates

perch_FA_covariates <- perch_FA_prop2[,c(1:3,51:59),]

perch_FA_totals <- perch_FA_prop2[,c(17,26,35,45)]

summary(perch_FA_totals)

# Pull out fatty acid composition matrix

perch_FA_prop_matrix <- perch_FA_prop2[,c(7:16,18:25,27:34,36:44)]

# Calculate distance matrix
perch_FA_diss <- as.matrix(vegdist(perch_FA_prop_matrix, 
                                   method = "euclidean", 
                                   na.rm = TRUE), 
                  labels = TRUE)

NMDS_scree(perch_FA_diss)  # 4 dimensions looks good

set.seed(5465)

perch_FA_nMDS1 <- 
  metaMDS(perch_FA_diss,
          distance = "euclidean",
          k = 4,
          trymax = 250,
          wascores = TRUE,
          expand = TRUE,
          autotransform = FALSE)

# Shepards test/goodness of fit
goodness(perch_FA_nMDS1)
stressplot(perch_FA_nMDS1)

perch_FA_data.scores <- as.data.frame(scores(perch_FA_nMDS1))
perch_FA_data.scores2 <- cbind(perch_FA_data.scores,
                               perch_FA_covariates)

perch_FA_scores <- `sppscores<-`(perch_FA_nMDS1, perch_FA_prop_matrix)

perch_FA_variable_scores <- 
  as.data.frame(perch_FA_scores$species)

perch_FA_variable_scores$variable <- 
  rownames(perch_FA_variable_scores)

## Generate hulls ----

perch_FA_data.scores2$MPconcentration <- 
  as.factor(perch_FA_data.scores2$MPconcentration)

perch_FA_hulls <- data.frame()

for(i in 1:length(unique(perch_FA_data.scores2$MPconcentration))) {
  hull <-
    perch_FA_data.scores2[perch_FA_data.scores2$MPconcentration ==
                            unique(perch_FA_data.scores2$MPconcentration)[i],
                          ][chull(perch_FA_data.scores2[perch_FA_data.scores2$MPconcentration ==
                                                          unique(perch_FA_data.scores2$MPconcentration)[i],
                                                        c(1:2)]),]
  perch_FA_hulls <- rbind(perch_FA_hulls, hull)
}

## Plot ----

png("Perch MDS Plot.png",
    width = 19,
    height= 15, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = perch_FA_hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   fill = MPconcentration,
                   colour = MPconcentration),
               alpha = 0.75,
               size = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_label(data = perch_FA_hulls,
             aes(x = NMDS1,
                 y = NMDS2,
                 label = ID),
             size = 4,
             alpha = 0.3) +
  geom_text(data = perch_FA_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = variable,
                angle = MDS1*MDS2*1.5),
            alpha = 0.9,
            size = 3,
            colour = "purple") +
  scale_fill_brewer(type = "seq",
                    palette = "YlGnBu",
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "seq",
                      palette = "YlGnBu",
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1

dev.off()

perch_FA_pca <- prcomp(perch_FA_prop_matrix, 
                      center = TRUE,
                      scale. = TRUE)

summary(perch_FA_pca)

library(ggbiplot)

ggbiplot(perch_FA_pca,
         labels = perch_FA_covariates$ID,
         ellipse = TRUE,
         groups = as.factor(perch_FA_covariates$MPconcentration))


# Essential fatty acids comparison ----

## Put data into long form ----

FA_prop_long <- 
  FAs_prop[,c(1:3,5,6,17,20,23,26,27,31,35,36,40,42,45:47,51)] %>%
  pivot_longer(names(FAs_prop[c(17,20,23,26,27,31,35,36,40,42,45:47)]),
               names_to = "metric",
               values_to = "value")

FA_prop_long$metric <- as.factor(FA_prop_long$metric)

FA_prop_long$sample.type <- 
  mapvalues(FA_prop_long$sample.type,
            from = levels(FA_prop_long$sample.type),
            to = c("Mysis Fish Food",
                   "Yellow Perch",
                   "Zooplankton"))

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
       y = "prop of Total Fatty Acids") +
  scale_colour_manual(values = c("orange", "purple", "blue"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

dev.off()

# Perch HMSC model----

## Specify data ----

perch_FA_prop_predictors <- 
  perch_FA_prop2[,c(51,54,64)]  # pull out predictors

perch_FA_prop_RE <- data.frame(corral = factor(perch_FA_prop2[,2]))
perch_FA_prop_rlevels <- HmscRandomLevel(units = perch_FA_prop_RE$corral)

perch_FA_prop_response <- 
  perch_FA_prop2[,c(7:16,18:25,27:34,36:44)]  # pull out FA proportions

## Specify model structure ----

model1 <- 
  Hmsc(Y = perch_FA_prop_response,  # response data
       XData = perch_FA_prop_predictors,  # covariates
       XFormula = ~ body.weight + MPconcentration + date2,  # model formula
       XScale = TRUE,  # scale covariates for fixed effects,
       studyDesign = perch_FA_prop_RE,
       ranLevels = list(corral = perch_FA_prop_rlevels),
       distr = "normal")

## Run MCMC chains ----

set.seed(6461)

perch_FA_prop_run1 <- sampleMcmc(model1,
                                    thin = 1,
                                    samples = 2000,
                                    transient = 100,
                                    nChains = 3,
                                    nParallel = 3,
                                    verbose = 1000)

## Check convergence ----

perch_FA_prop_post1 <- convertToCodaObject(perch_FA_prop_run1)

effectiveSize(perch_FA_prop_post1$Beta)
gelman.diag(perch_FA_prop_post1$Beta, 
            transform = TRUE,
            multivariate = FALSE)$psrf

plot(perch_FA_prop_post1$Beta)


## Assess explanatory power ----

perch_FA_prop_pred1 <- computePredictedValues(perch_FA_prop_run1)
evaluateModelFit(hM = perch_FA_prop_run1, predY = perch_FA_prop_pred1)


## Cross validation ----

partition1 <- createPartition(perch_FA_prop_run1, nfolds = 2)
perch_FA_prop_pred1.1 <- 
  computePredictedValues(perch_FA_prop_run1, partition = partition1)

evaluateModelFit(hM = perch_FA_prop_run1, 
                 predY = perch_FA_prop_pred1.1)

## Look at slope estimates ----

perch_FA_prop_postBeta <- 
  getPostEstimate(perch_FA_prop_run1, parName = "Beta")
plotBeta(perch_FA_prop_run1, post = perch_FA_prop_postBeta, 
         param = "Support", supportLevel = 0.95)

# Zooplankton HMSC model----

## Specify data ----

zoop_FA_prop$datefactor <- as.factor(zoop_FA_prop$date)

zoop_FA_prop_predictors <- zoop_FA_prop[,c(51,53)]  # pull out predictor

zoop_FA_prop_RE <- data.frame(corral = zoop_FA_prop[,2])
zoop_FA_prop_rlevels <- HmscRandomLevel(units = zoop_FA_prop_RE$corral)

zoop_FA_prop_response <- 
  zoop_FA_prop[,c(7:16,18:25,27:34,36:44)]  # pull out FA proportions

## Specify model structure ----

zoop_FA_prop_model1 <- 
  Hmsc(Y = zoop_FA_prop_response,  # response data
       XData = zoop_FA_prop_predictors,  # covariates
       XFormula = ~ MPconcentration * datefactor,  # model formula
       XScale = TRUE,  # scale covariates for fixed effects,
       studyDesign = zoop_FA_prop_RE,
       ranLevels = list(corral = zoop_FA_prop_rlevels),
       distr = "normal")

## Run MCMC chains ----

set.seed(6561)

zoop_FA_prop_run1 <- sampleMcmc(zoop_FA_prop_model1,
                                    thin = 1,
                                    samples = 2000,
                                    transient = 100,
                                    nChains = 3,
                                    nParallel = 3,
                                    verbose = 1000)

## Check convergence ----

zoop_FA_prop_post1 <- convertToCodaObject(zoop_FA_prop_run1)

effectiveSize(zoop_FA_prop_post1$Beta)
gelman.diag(zoop_FA_prop_post1$Beta, 
            transform = TRUE,
            multivariate = FALSE)$psrf

plot(zoop_FA_prop_post1$Beta)


## Assess explanatory power ----

zoop_FA_prop_pred1 <- computePredictedValues(zoop_FA_prop_run1)
evaluateModelFit(hM = zoop_FA_prop_run1, predY = zoop_FA_prop_pred1)


## Cross validation ----

partition1 <- createPartition(zoop_FA_prop_run1, nfolds = 2)
zoop_FA_prop_pred1.1 <- 
  computePredictedValues(zoop_FA_prop_run1, partition = partition1)

evaluateModelFit(hM = zoop_FA_prop_run1, 
                 predY = zoop_FA_prop_pred1.1)

## Look at slope estimates ----

zoop_FA_prop_postBeta <- 
  getPostEstimate(zoop_FA_prop_run1, parName = "Beta")
plotBeta(zoop_FA_prop_run1, post = zoop_FA_prop_postBeta, 
         param = "Support", supportLevel = 0.95)

