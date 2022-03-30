# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(plyr)
library(dplyr)
library(MuMIn)
library(vegan)
library(tidyr)
library(DirichletReg)

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



# Analyses ----

## Perch analyses ----

### nMDS ----

# Pull out covariates

perch_FA_covariates <- perch_FA_prop2[,c(1:3,51:66),]

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

#### Generate hulls ----

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

#### Plot ----

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
  scale_fill_brewer(type = "div",
                    palette = "RdYlGn",
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1

dev.off()


### Dirichlet regression ----

# Re-scale response so it sums to 1

perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] 

adjusted_Y <- 
  DR_data(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] / 
            rowSums(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)]))

# Run model

perch_FA_dir_mod1 <-
  DirichReg(adjusted_Y ~ scaled.body.weight + log(MPconcentration + 1)|
              corral,
            model = "alternative",
            data = perch_FA_prop2)

perch_FA_dir_mod1
summary(perch_FA_dir_mod1)

# Predict from original data

perch_FA_dir_predict1 <- predict(perch_FA_dir_mod1, mu = TRUE, phi = TRUE)

perch_FA_prop2$precision <- perch_FA_dir_predict1$phi

#### Plot precision ----

ggplot(perch_FA_prop2) +
  geom_point(aes(x = MPconcentration,
                 y = precision,
                 colour = corral),
             size = 4) +
labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
     y = "Precision") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_colour_brewer(type = "qual",
                    palette = 1,
                    name = "Corral") +
  theme1

# Predict from new data

perch_FA_dir_newdata <- perch_FA_prop2

perch_FA_dir_newdata$scaled.body.weight = 0

perch_FA_dir_newdata$scaled.date = 0

perch_FA_dir_newdata$scaled.pop = 0

perch_FA_dir_newpredict <- predict(perch_FA_dir_mod1, 
                                   newdata = perch_FA_dir_newdata,
                                   mu = TRUE, phi = FALSE)

colnames(perch_FA_dir_newpredict) <- 
  colnames(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)])

perch_FA_dir_newdata$MPconcentration <-
  (perch_FA_dir_newdata$scaled.MPconcentration * 
     sd(perch_FA_prop2$MPconcentration)) +
  mean(perch_FA_prop2$MPconcentration)

# Put predictions into long form

perch_FA_dir_newpredict <-
  cbind(perch_FA_dir_newdata[,c(1:6,51:68)], perch_FA_dir_newpredict)

perch_FA_dir_newpredict_long <- 
  perch_FA_dir_newpredict %>%
    pivot_longer(names(perch_FA_dir_newpredict[c(25:59)]),
               names_to = "metric",
               values_to = "value")

perch_FA_dir_newpredict_long$metric <- 
  as.factor(perch_FA_dir_newpredict_long$metric)

# Put original data into long form

perch_FA_prop_long <- 
  perch_FA_prop2[,c(1:3,5,6:16,18:25,27:34,36:44,51:53)] %>%
  pivot_longer(names(zoop_FA_prop[c(7:16,18:25,27:34,36:44)]),
               names_to = "metric",
               values_to = "value")

perch_FA_prop_long$metric <- as.factor(perch_FA_prop_long$metric)


#### Plot predictions ----

ggplot() +
  geom_line(data = perch_FA_dir_newpredict_long,
            aes(x = MPconcentration,
                y = value),
            colour = "red") +
  geom_point(data = perch_FA_prop_long,
             aes(x = MPconcentration,
                 y = value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion Fatty Acid") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ metric, ncol = 5, scales = "free_y") +
  theme1



## Zooplankton analyses ----

### nMDS ----

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

#### Generate hulls ----

zoop_FA_data.scores2$MPconcentration <- 
  as.factor(zoop_FA_data.scores2$MPconcentration)

zoop_FA_data.scores2$date <- 
  as.factor(zoop_FA_data.scores2$date)

zoop_FA_hulls <- data.frame()

for(i in 1:length(unique(zoop_FA_data.scores2$date))) {
  hull <-
    zoop_FA_data.scores2[zoop_FA_data.scores2$date ==
                           unique(zoop_FA_data.scores2$date)[i],
    ][chull(zoop_FA_data.scores2[zoop_FA_data.scores2$date ==
                                   unique(zoop_FA_data.scores2$date)[i],
                                 c(1:2)]),]
  zoop_FA_hulls <- rbind(zoop_FA_hulls, hull)
}

### Plot ----

png("Zooplankton MDS Plot.png",
    width = 19,
    height= 13, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = zoop_FA_hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   colour = date),
               alpha = 0.05,
               size = 1,
               fill = "black") +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = zoop_FA_hulls,
             aes(x = NMDS1,
                 y = NMDS2,
                 fill = MPconcentration),
             alpha = 0.75,
             size = 4,
             shape = 21) +
  geom_text(data = zoop_FA_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = variable,
                angle = MDS1*MDS2*5000),
            alpha = 0.9,
            size = 2,
            colour = "purple") +
  scale_colour_brewer(type = "qual",
                    palette = "Set1",
                    direction = -1,
                    name = "Date") +
  scale_fill_brewer(type = "div",
                      palette = "RdYlGn",
                    direction = -1,
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1

dev.off()

### Dirichlet regression ----

# Subset endpoint data

zoop_FA_prop_end <- subset(zoop_FA_prop, date > "2021-08-01" )

# Re-scale response so it sums to 1

zoop_FA_prop_end[,c(7:16,18:25,27:34,36:44)] 

adjusted_Y <- 
  DR_data(zoop_FA_prop_end[,c(7:16,18:25,27:34,36:44)] / 
            rowSums(zoop_FA_prop_end[,c(7:16,18:25,27:34,36:44)]))

# Run model

zoop_FA_dir_mod1 <-
  DirichReg(adjusted_Y ~ log(MPconcentration + 1)|
              log(MPconcentration + 1),
            model = "alternative",
            data = zoop_FA_prop_end)

zoop_FA_dir_mod1
summary(zoop_FA_dir_mod1)

# Predict from original data

zoop_FA_dir_predict1 <- predict(zoop_FA_dir_mod1, mu = TRUE, phi = FALSE)

zoop_FA_prop_end_predict <- 
  cbind(zoop_FA_prop_end[,c(1:6,51:54)],
        zoop_FA_dir_predict1)

# Put predictions into long form

zoop_FA_prop_end_predict_long <- 
  zoop_FA_prop_end_predict %>%
  pivot_longer(names(zoop_FA_prop_end_predict[c(11:45)]),
               names_to = "metric",
               values_to = "value")

zoop_FA_prop_end_predict_long$metric <- 
  as.factor(zoop_FA_prop_end_predict_long$metric)

# Put original data into long form

zoop_FA_prop_long <- 
  zoop_FA_prop[,c(1:3,5,6:16,18:25,27:34,36:44,51:54)] %>%
  pivot_longer(names(zoop_FA_prop[c(7:16,18:25,27:34,36:44)]),
               names_to = "metric",
               values_to = "value")

zoop_FA_prop_long$metric <- as.factor(zoop_FA_prop_long$metric)

zoop_FA_prop_long_end <- 
  subset(zoop_FA_prop_long, date > "2021-08-01" )

#### Plot predictions ----

ggplot() +
  geom_line(data = zoop_FA_prop_end_predict_long,
            aes(x = MPconcentration,
                y = value),
            colour = "red") +
  geom_point(data = zoop_FA_prop_long_end,
             aes(x = MPconcentration,
                 y = value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion Fatty Acid") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ metric, ncol = 5, scales = "free_y") +
  theme1

# DHA prediction

zoop_FA_prop_end$DHA.pred <- zoop_FA_dir_predict1[,33]

ggplot(zoop_FA_prop_end) +
  geom_line(aes(x = MPconcentration,
                y = DHA.pred)) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion DHA") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

# EPA prediction

zoop_FA_prop_end$EPA.pred <- zoop_FA_dir_predict1[,31]

ggplot(zoop_FA_prop_end) +
  geom_line(aes(x = MPconcentration,
                y = EPA.pred)) +
  geom_point(aes(x = MPconcentration,
                 y = C_22.6n.3)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion EPA") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1





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




perch_FA_pca <- prcomp(perch_FA_prop_matrix, 
                      center = TRUE,
                      scale. = TRUE)

summary(perch_FA_pca)

library(ggbiplot)

ggbiplot(perch_FA_pca,
         labels = perch_FA_covariates$ID,
         ellipse = TRUE,
         groups = as.factor(perch_FA_covariates$MPconcentration))


