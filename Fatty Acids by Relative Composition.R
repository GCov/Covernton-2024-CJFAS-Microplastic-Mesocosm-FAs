# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)
library(vegan)
library(tidyr)

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

# Separate out by sample type ----

perch_FA_percent <- subset(FAs_percent, sample.type == "perch")
zoop_FA_percent <- subset(FAs_percent, sample.type == "zooplankton")
food_FA_percent <- subset(FAs_percent, sample.type == "fish food")

# Assign corral a MP concentration of 0 for zoop_FA

zoop_FA_percent$MPconcentration[zoop_FA_percent$corral == "A"] <- 
  24

# Add in perch biometrics data ----

# Combine with perch FA data

perch_FA_percent2 <- left_join(perch_FA_percent, 
                               perch_biometrics, 
                       by = c("ID", "corral"))

## Plot by fatty acid composition ----

# Put data into long form and remove totals

perch_FA_percent_long <- 
  perch_FA_percent2[,c(1:3,6,7:16,18:25,27:34,36:44,51:59)] %>%
  pivot_longer(names(perch_FA_percent2[c(7:16,18:25,27:34,36:44)]),
               names_to = "fatty.acid",
               values_to = "percent")

# Plot

ggplot(perch_FA_percent_long) +
  geom_col(aes(x = ID,
               y = percent,
               fill = fatty.acid),
           colour = "black",
           size = 0.25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = expression(paste("MP Exposure Concentration (particles "~L^-1*")")),
       y = "Percent Composition") +
  scale_fill_hue(name = "Fatty Acid",
                 l = 75) +
  facet_grid(.~MPconcentration,
             scales = "free_x",
             switch = "x",
             space = "free_x") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Zooplankton nMDS ----

# Pull out covariates

zoop_FA_covariates <- zoop_FA_percent[,c(1:3,51),]

zoop_FA_totals <- zoop_FA_percent[,c(17,26,35,45)]

# Pull out fatty acid composition matrix

zoop_FA_percent_matrix <- zoop_FA_percent[,c(7:16,18:25,27:34,36:44)]

# Calculate distance matrix
zoop_FA_diss <- as.matrix(vegdist(zoop_FA_percent_matrix, 
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

zoop_FA_scores <- `sppscores<-`(zoop_FA_nMDS1, zoop_FA_percent_matrix)

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

zoop_FA_pca <- prcomp(zoop_FA_percent_matrix, 
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

perch_FA_covariates <- perch_FA_percent2[,c(1:3,51:59),]

perch_FA_totals <- perch_FA_percent2[,c(17,26,35,45)]

summary(perch_FA_totals)

# Pull out fatty acid composition matrix

perch_FA_percent_matrix <- perch_FA_percent2[,c(7:16,18:25,27:34,36:44)]

# Calculate distance matrix
perch_FA_diss <- as.matrix(vegdist(perch_FA_percent_matrix, 
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

perch_FA_scores <- `sppscores<-`(perch_FA_nMDS1, perch_FA_percent_matrix)

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

perch_FA_pca <- prcomp(perch_FA_percent_matrix, 
                      center = TRUE,
                      scale. = TRUE)

summary(perch_FA_pca)

library(ggbiplot)

ggbiplot(perch_FA_pca,
         labels = perch_FA_covariates$ID,
         ellipse = TRUE,
         groups = as.factor(perch_FA_covariates$MPconcentration))
