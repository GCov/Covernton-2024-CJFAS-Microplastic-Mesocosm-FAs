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
library(easyCODA)

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
       y = "Proportion of Total Fatty Acids") +
  scale_colour_manual(values = c("orange", "purple", "blue"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  theme1

dev.off()



# Analyses ----

## Perch analyses ----

### PCA ----

#### Pull out covariates ----

perch_FA_prop_covariates <- perch_FA_prop2[,c(1:3,51:67),]

perch_FA_prop_totals <- perch_FA_prop2[,c(17,26,35,45)]

summary(perch_FA_prop_totals)

#### Pull out fatty acid composition matrix ----

perch_FA_prop_matrix <- perch_FA_prop2[,c(7:16,18:25,27:34,36:44)]

#### Run PCA ----
perch_FA_prop_pca <- rda(perch_FA_prop_matrix,
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
perch_FA_prop_PCA_species$FA <- FA.names

perch_FA_prop_PCA_site <- cbind(perch_FA_prop_covariates,
                                perch_FA_prop_PCA_site[, 1:2])

#### Plot ----

png("Perch FA Proportions PCA.png",
    width = 19,
    height= 14, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_segment(data = perch_FA_prop_PCA_species,
               aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_point(data = perch_FA_prop_PCA_site,
             aes(x = PC1,
                 y = PC2,
                 colour = as.factor(MPconcentration)),
             size = 4,
             alpha = 0.75) +
  geom_text(data = perch_FA_prop_PCA_species,
            aes(x = PC1, 
                y = PC2, 
                label = FA),
            alpha = 0.5,
            size = 3,
            colour = "purple") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      name = 
                        expression(paste("Exposure propentration (MPs"~L^-1*")")),
                      direction = -1) +
  labs(x = "PC1",
       y = "PC2") +
  theme1

dev.off()


### nMDS ----

# Remove any FAs that contribute less than 0.01% and re-scale response so it 
# sums to 1

trimmed_perch_FA <- 
  data.frame(perch_FA_prop2[,c(7:16,18:25,27:34,36:44)] %>% 
               select(where(function(x){mean(x) >= 0.01})))

trimmed_perch_FA <- 
  trimmed_perch_FA / rowSums(trimmed_perch_FA)

# Calculate distance matrix
perch_FA_prop_diss <- as.matrix(vegdist(trimmed_perch_FA, 
                                   method = "bray", 
                                   na.rm = TRUE), 
                           labels = TRUE)

NMDS_scree(perch_FA_prop_diss)  # 4 dimensions looks good

set.seed(5465)

perch_FA_prop_nMDS1 <- 
  metaMDS(perch_FA_prop_diss,
          distance = "bray",
          k = 4,
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
  FA.names[c(3,5,7,13,15,16,19,23,26:28,31:33)]

#### Generate hulls ----

perch_FA_prop_data.scores2$MPconcentration <- 
  as.factor(perch_FA_prop_data.scores2$MPconcentration)

perch_FA_prop_hulls <- data.frame()

for(i in 1:length(unique(perch_FA_prop_data.scores2$MPconcentration))) {
  hull <-
    perch_FA_data.scores2[perch_FA_prop_data.scores2$MPconcentration ==
                            unique(perch_FA_prop_data.scores2$MPconcentration)[i],
    ][chull(perch_FA_data.scores2[perch_FA_prop_data.scores2$MPconcentration ==
                                    unique(perch_FA_prop_data.scores2$MPconcentration)[i],
                                  c(1:2)]),]
  perch_FA_prop_hulls <- rbind(perch_FA_prop_hulls, hull)
}

#### Plot ----

png("Perch FA Proportions nMDS Plot.png",
    width = 19,
    height= 14, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = perch_FA_prop_hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   fill = MPconcentration,
                   colour = MPconcentration),
               alpha = 0.3,
               size = 0.5) +
  geom_segment(data = perch_FA_prop_variable_scores,
               aes(x = 0, y = 0, xend = MDS1*0.9, yend = MDS2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = perch_FA_prop_hulls,
             aes(x = NMDS1,
                 y = NMDS2,
                 colour = MPconcentration),
             size = 4,
             alpha = 0.5) +
  geom_text(data = perch_FA_prop_variable_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = FA),
            alpha = 0.9,
            size = 3,
            colour = "purple") +
  scale_fill_brewer(type = "div",
                    palette = "RdYlGn",
                    direction = -1,
                    name = 
                      expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      direction = -1,
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  theme1

dev.off()

### PERMANOVA ####

perch.permanova <- 
  adonis2(perch_FA_prop_diss ~ log(MPconcentration + 1) + body.weight,
          permutations = 999,
          data = perch_FA_prop_covariates,
          method = "bray",
          parallel = 8,
          strata = perch_FA_prop_covariates$corral)

perch.permanova  

# No significant effect of MP concentration or body weight


### Dirichlet regression ----

DR_Y <- DR_data(trimmed_FA)

# Run model

perch_FA_dir_mod1 <-
  DirichReg(DR_Y ~ body.weight + log(MPconcentration + 1)|
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

perch_FA_dir_newdata$body.weight = mean(perch_FA_prop2$body.weight)

perch_FA_dir_newpredict <- predict(perch_FA_dir_mod1, 
                                   newdata = perch_FA_dir_newdata,
                                   mu = TRUE, phi = FALSE)

colnames(perch_FA_dir_newpredict) <- 
  colnames(trimmed_FA)

perch_FA_dir_newdata$MPconcentration <-
  (perch_FA_dir_newdata$scaled.MPconcentration * 
     sd(perch_FA_prop2$MPconcentration)) +
  mean(perch_FA_prop2$MPconcentration)

# Put predictions into long form

perch_FA_dir_newpredict <-
  cbind(perch_FA_dir_newdata[,c(1:6,51:68)], perch_FA_dir_newpredict)

perch_FA_dir_newpredict_long <- 
  perch_FA_dir_newpredict %>%
    pivot_longer(names(perch_FA_dir_newpredict[c(25:39)]),
               names_to = "metric",
               values_to = "value")

perch_FA_dir_newpredict_long$metric <- 
  as.factor(perch_FA_dir_newpredict_long$metric)

# Put original data into long form

perch_FA_prop_long <- 
  trimmed_FA %>%
  pivot_longer(names(trimmed_FA),
               names_to = "metric",
               values_to = "value")

perch_FA_prop_long$metric <- as.factor(perch_FA_prop_long$metric)


#### Plot predictions ----

png("Perch Dirichlet Plot.png",
    width = 29,
    height= 15, 
    units = "cm",
    res = 600)

ggplot() +
  geom_line(data = perch_FA_dir_newpredict_long,
            aes(x = MPconcentration,
                y = value),
            colour = "red") +
  geom_point(data = perch_FA_prop_long,
             aes(x = perch_FA_dir_newpredict_long$MPconcentration,
                 y = value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion Fatty Acid") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ metric, ncol = 5) +
  theme1

dev.off()

### Explore different indicators ----

#### PUFAs ----

perch_FA_prop_PUFA_mod1 <- 
  glmmTMB(PUFAs ~ scaled.MPconcentration + scaled.body.weight,
          data = perch_FA_prop2,
          family = beta_family(link = "logit"))

plot(simulateResiduals(perch_FA_prop_PUFA_mod1))

summary(perch_FA_prop_PUFA_mod1)

#### SFAs ----

perch_FA_prop_SFA_mod1 <- 
  glmmTMB(total_SFAs ~ scaled.MPconcentration,
          data = perch_FA_prop2,
          family = beta_family(link = "logit"))

plot(simulateResiduals(perch_FA_prop_SFA_mod1))

summary(perch_FA_prop_SFA_mod1)

### Try iterative process from Graeve and Greenacre (2020) ----

# Average proportions
round(colMeans(perch_FA_prop_matrix), 3)

### Normalize the data
perch_FA_prop_matrix2 <- perch_FA_prop_matrix / rowSums(perch_FA_prop_matrix)

### Define function for output of ratios
print.ratios <- function(rationames, R2, procr=NA, N=10) {
  # function prints ratios and the corresponding R2, optionally Procrustes correlations
  # prints first 10 ratios by default  
  # split names of parts in ratio into FA1 and FA2
  # notice that the ratios can be reported as FA1/FA2 or FA2/FA1  
  foo    <- strsplit(rationames, "/")
  parts  <- matrix(unlist(foo), ncol=2, byrow=TRUE)
  df   <- as.data.frame(parts)[1:N,]
  if(is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2))
    colnames(df) <- c("FA1", "FA2","R2")
  }
  if(!is.na(procr)) {
    df <- cbind(df, round(100*R2[1:N], 2), round(procr[1:N], 3))
    colnames(df) <- c("FA1", "FA2", "R2","Procrustes")
  }
  print(df[1:N,])
}

#### Step 1 ----
FA.step1 <- STEP(perch_FA_prop_matrix2, nsteps = 1, top=20)
print.ratios(FA.step1$names.top, FA.step1$R2.top)

# Select 1st option: 16:1n-7/22:6n-3 (69.6% variance)

# logratios1, 2, etc... will gather the logratio values chosen
# numratios1, 2, etc... will gather the numbers of the two parts 
logratios1 <- FA.step1$logratios.top[,1]
numratios1 <- FA.step1$ratios.top[1,]

#### Step 2 ----
FA.step2 <- STEP(perch_FA_prop_matrix2, nsteps = 1, top=20, 
                 previous=logratios1)
print.ratios(FA.step2$names.top, FA.step2$R2.top)

# Select 1st option: 14:0/18:3n-3 (82.3% variance)

### update logratios1 to logratios2, and numratios1 to numratios2
logratios2 <- cbind(logratios1, FA.step2$logratios.top[,1])
numratios2 <- rbind(numratios1, FA.step2$ratios.top[1,])

#### Step 3 ----

FA.step3 <- STEP(perch_FA_prop_matrix2, nsteps = 1, top=20, previous=logratios2)
print.ratios(FA.step3$names.top, FA.step3$R2.top)

# Select 3rd option: 16:1n-7/20:4n-6 (85.7% variance)

logratios3 <- cbind(logratios2, FA.step3$logratios.top[,3])
numratios3 <- rbind(numratios2, FA.step3$ratios.top[3,])

#### Step 4 ----

FA.step4 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios3)
print.ratios(FA.step4$names.top, FA.step4$R2.top)

# Select 4th option: 16:1n-7/18:4n-3 (88.5% variance)

logratios4 <- cbind(logratios3, FA.step4$logratios.top[,4])
numratios4 <- rbind(numratios3, FA.step4$ratios.top[4,])

#### Step 5 ----
FA.step5 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios4)
print.ratios(FA.step5$names.top, FA.step5$R2.top)

# Select 9th option: 16:1n-7/20:5n-3 (91.1% variance)

logratios5 <- cbind(logratios4, FA.step5$logratios.top[,9])
numratios5 <- rbind(numratios4, FA.step5$ratios.top[9,])

#### Step 6 ----
FA.step6 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios5)
print.ratios(FA.step6$names.top, FA.step6$R2.top, N = 20)

# Select 2nd option: 18:2n-6/22:6n-3 (93.4% variance)

logratios6 <- cbind(logratios5, FA.step6$logratios.top[,2])
numratios6 <- rbind(numratios5, FA.step6$ratios.top[2,])

#### Step 7 ----
FA.step7 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios6)
print.ratios(FA.step7$names.top, FA.step7$R2.top, N = 20)

# Select 1st option: 16:0/18:2n-6 (95.0% variance)

logratios7 <- cbind(logratios6, FA.step7$logratios.top[,1])
numratios7 <- rbind(numratios6, FA.step7$ratios.top[1,])

#### Step 8 ----
FA.step8 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios7)
print.ratios(FA.step8$names.top, FA.step8$R2.top, N = 20)

# Select 1st option: 18:2n-6/18:4n-3 (100.0% variance)

logratios8 <- cbind(logratios7, FA.step8$logratios.top[,1])
numratios8 <- rbind(numratios7, FA.step8$ratios.top[1,])

#### Step 9 ----
FA.step9 <- STEP(perch_FA_prop_matrix2, nsteps=1, top=20, previous=logratios8)
print.ratios(FA.step9$names.top, FA.step9$R2.top, N = 20)

## DONE - Stop at Step 8

### These are the ratios chosen in the 8 steps (numbers first, then names)
rownames(numratios8) <- paste("Step",1:8,sep="")
colnames(numratios8) <- c("FA1","FA2")
finalratios <- 
  as.data.frame(cbind(numratios8, 
                      Ratio=paste(colnames(perch_FA_prop_matrix2)[numratios8[, 1]],
                                  "/",colnames(perch_FA_prop_matrix2)[numratios8[,2]],
                                  sep="")))
finalratios

colnames(logratios8) <- finalratios[,3]

# The 9 parts used in the 8 ratios

partsinratios <- sort(unique(as.numeric(numratios8)))
colnames(perch_FA_prop_matrix2)[partsinratios]

### the LRA of the full data set 
### low-contributing FAs are de-accentuated by plotting them in light red
### from the LRA object
FA0.lra <- LRA(perch_FA_prop_matrix2)
### .ccc are the column contribution coordinates
### .ctr TRUE for high contributor, otherwise FALSE
FA0.ccc <- FA0.lra$colcoord * sqrt(FA0.lra$colmass)
FA0.ctr <- (FA0.ccc[,1]^2 * FA0.lra$sv[1]^2 + FA0.ccc[,2]^2 * FA0.lra$sv[2]^2) / 
  (FA0.lra$sv[1]^2 + FA0.lra$sv[2]^2) > 1/ncol(perch_FA_prop_matrix2)
### only show parts for which FA0.ctr = TRUE (high contributors)
FA0.lra$colcoord <- FA0.lra$colcoord[FA0.ctr,]
FA0.lra$colmass  <- FA0.lra$colmass[FA0.ctr]
FA0.lra$colnames <- FA0.lra$colnames[FA0.ctr]
PLOT.LRA(FA0.lra, map="contribution")
text(FA0.ccc[!FA0.ctr,], 
     labels=colnames(perch_FA_prop_matrix2)[!FA0.ctr], 
     col="pink", cex=0.6)

#### Plot PCA ----

### Run PCA ----
perch_FA_reduced_pca <- rda(logratios8,
                            scale. = FALSE)

# Bar plot of relative eigenvalues
barplot(as.vector(perch_FA_reduced_pca$CA$eig)/sum(perch_FA_reduced_pca$CA$eig))

# Calculate percentage of variance explained by first 2 aaxes
sum((as.vector(perch_FA_reduced_pca$CA$eig)/sum(perch_FA_reduced_pca$CA$eig))[1:2])

# 'Site' scores
perch_FA_reduced_pca_site <- perch_FA_reduced_pca$CA$u

# 'Species' scores
perch_FA_reduced_pca_species <- data.frame(perch_FA_reduced_pca$CA$v)
perch_FA_reduced_pca_species$LR <- rownames(perch_FA_reduced_pca_species)

perch_FA_reduced_pca_site <- cbind(perch_FA_prop_covariates,
                                   perch_FA_reduced_pca_site[, 1:2])

perch_FA_reduced_envfit <- 
  envfit(perch_FA_reduced_pca, 
         perch_FA_prop_covariates[,c(4, 7)])

perch_FA_reduced_predictors <- 
  as.data.frame(perch_FA_reduced_envfit$vectors$arrows)

perch_FA_reduced_predictors$labels <- c("MP concentration",
                                        "Body weight")

### Plot ----

png("Perch FA Log-ratio PCA.png",
    width = 19,
    height= 14, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_segment(data = perch_FA_reduced_pca_species,
               aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_segment(data = perch_FA_reduced_predictors,
               aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "blue",
               alpha = 0.5) +
  geom_point(data = perch_FA_reduced_pca_site,
             aes(x = PC1,
                 y = PC2,
                 colour = as.factor(MPconcentration)),
             size = 4,
             alpha = 0.75) +
  geom_text(data = perch_FA_reduced_pca_species,
            aes(x = PC1, 
                y = PC2, 
                label = LR),
            alpha = 0.5,
            size = 3,
            colour = "purple") +
  geom_text(data = perch_FA_reduced_predictors,
            aes(x = PC1, 
                y = PC2, 
                label = labels),
            alpha = 0.5,
            size = 3,
            colour = "blue") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      name = 
                        expression(paste("Exposure propentration (MPs"~L^-1*")")),
                      direction = -1) +
  labs(x = "PC1",
       y = "PC2") +
  theme1

dev.off()


## Zooplankton analyses ----

### PCA ----

#### Pull out covariates ----

zoop_FA_prop_covariates <- zoop_FA_prop[,c(1:3,51:54),]

zoop_FA_prop_totals <- zoop_FA_prop[,c(17,26,35,45)]

summary(zoop_FA_prop_totals)

# Pull out fatty acid composition matrix

zoop_FA_prop_matrix <- as.matrix(zoop_FA_prop[,c(7:16,18:25,27:34,36:44)])

#### Run PCA ----
zoop_FA_prop_pca <- rda(zoop_FA_prop_matrix,
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
zoop_FA_prop_PCA_species$FA <- FA.names

zoop_FA_prop_PCA_site <- cbind(zoop_FA_prop_covariates,
                               zoop_FA_prop_PCA_site[, 1:2])

#### Plot ----

png("Zooplankton FA Proportions PCA.png",
    width = 19,
    height= 14, 
    units = "cm",
    res = 600)

ggplot() +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_segment(data = zoop_FA_prop_PCA_species,
               aes(x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_point(data = zoop_FA_prop_PCA_site,
             aes(x = PC1,
                 y = PC2,
                 colour = as.factor(MPconcentration),
                 shape = as.factor(date)),
             size = 4,
             alpha = 0.75) +
  geom_text(data = zoop_FA_prop_PCA_species,
            aes(x = PC1, 
                y = PC2, 
                label = FA),
            alpha = 0.5,
            size = 3,
            colour = "purple") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      name = 
                        expression(paste("Exposure propentration (MPs"~L^-1*")")),
                      direction = -1) +
  scale_shape(name = "Date") +
  labs(x = "PC1",
       y = "PC2") +
  theme1

dev.off()

### nMDS ----

# Remove any FAs that contribute less than 0.01% and re-scale response so it 
# sums to 1

trimmed_zoop_FA <- 
  data.frame(zoop_FA_prop[,c(7:16,18:25,27:34,36:44)]) %>% 
  select(where(function(x){mean(x) >= 0.01}))

trimmed_zoop_FA <- trimmed_zoop_FA/rowSums(trimmed_zoop_FA)

# Calculate distance matrix
zoop_FA_prop_diss <- as.matrix(vegdist(trimmed_zoop_FA, 
                                  method = "bray", 
                                  na.rm = TRUE), 
                          labels = TRUE)

NMDS_scree(zoop_FA_prop_diss)  # 4 dimensions looks good

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

zoop_FA_variable_prop_scores$FA <- FA.names[c(3,5,7,13,15,16,18,19,20,23,26:28,31,33)]

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
    res = 600)

ggplot() +
  geom_polygon(data = zoop_FA_prop_hulls,
               aes(x = MDS1,
                   y = MDS2,
                   fill = date),
               colour = "black",
               alpha = 0.3,
               size = 0.5) +
  geom_segment(data = zoop_FA_prop_variable_scores,
               aes(x = 0, y = 0, xend = MDS1*0.9, yend = MDS2*0.9),
               arrow = arrow(type = "closed",
                             length = unit(0.2, "cm")),
               colour = "purple",
               alpha = 0.5) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed") +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(data = zoop_FA_prop_hulls,
             aes(x = MDS1,
                 y = MDS2,
                 shape = date,
                 colour = MPconcentration),
             alpha = 0.75,
             size = 3) +
  geom_text(data = zoop_FA_variable_prop_scores,
            aes(x = MDS1, 
                y = MDS2, 
                label = FA),
            alpha = 0.9,
            size = 8 / .pt,
            colour = "purple") +
  scale_colour_brewer(type = "div",
                      palette = "RdYlGn",
                      direction = -1,
                      name = 
                        expression(paste("Exposure Concentration (MPs"~L^-1*")"))) +
  scale_fill_brewer(type = "qual",
                    palette = "Set2",
                    direction = -1,
                      name = "Date") +
  scale_shape(name = "Date") +
  theme1

dev.off()

### PERMANOVA ####

zoop.permanova <- 
  adonis2(zoop_FA_prop_diss ~ log(MPconcentration + 1) * scaled.date,
          permutations = 999,
          data = zoop_FA_prop_covariates,
          method = "bray",
          parallel = 8,
          strata = zoop_FA_prop_covariates$corral)

zoop.permanova  

### Dirichlet regression ----

# Re-scale response so it sums to 1

trimmed_zoop_FA <- 
  data.frame(zoop_FA_prop[,c(7:16,18:25,27:34,36:44)]) %>% 
  select(where(function(x){mean(x) >= 0.01}))

trimmed_zoop_FA <- trimmed_zoop_FA/rowSums(trimmed_zoop_FA)

DR_zoop_Y <- DR_data(trimmed_zoop_FA)

# Run model

zoop_FA_dir_mod1 <-
  DirichReg(DR_zoop_Y ~ log(MPconcentration + 1) * scaled.date|
              scaled.date,
            model = "alternative",
            data = zoop_FA_prop)

zoop_FA_dir_mod1
summary(zoop_FA_dir_mod1)

# Predict from original data

zoop_FA_dir_predict1 <- 
  exp(predict(zoop_FA_dir_mod1, mu = TRUE, phi = FALSE)) - 1

zoop_FA_prop_end_predict <- 
  cbind(zoop_FA_prop[,c(1:6,51:54)],
        zoop_FA_dir_predict1)

# Put predictions into long form

zoop_FA_prop_end_predict_long <- 
  zoop_FA_prop_end_predict %>%
  pivot_longer(names(trimmed_zoop_FA),
               names_to = "metric",
               values_to = "value")

zoop_FA_prop_end_predict_long$metric <- 
  as.factor(zoop_FA_prop_end_predict_long$metric)

# Put original data into long form

zoop_FA_prop_long <- 
  cbind(zoop_FA_prop[,c(1:6,51:54)],
        trimmed_zoop_FA) %>%
  pivot_longer(names(trimmed_zoop_FA),
               names_to = "metric",
               values_to = "value")


#### Plot predictions ----

png("Zooplankton Dirichlet Plot.png",
    width = 19,
    height= 15, 
    units = "cm",
    res = 600)

ggplot() +
  geom_line(data = zoop_FA_prop_end_predict_long,
            aes(x = MPconcentration,
                y = value),
            colour = "red") +
  geom_point(data = zoop_FA_prop_long,
             aes(x = zoop_FA_prop_end_predict_long$MPconcentration,
                 y = value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Proportion Fatty Acid") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_grid(date ~ metric) +
  theme1

dev.off()

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


# Fish food ----

png("Fish Food Essential FAs Proportionss Plot.png",
    width = 9,
    height= 9, 
    units = "cm",
    res = 600)

ggplot(food_FA_prop) +
  geom_violin(aes(x = 1,
                  y = C_18.2n.6,
                  fill = "LA"),
              size = 0.25) +
  geom_violin(aes(x = 2,
                  y = C_20.4n.6,
                  fill = "ARA"),
              size = 0.25) +
  geom_violin(aes(x = 3,
                  y = C_18.3n.3,
                  fill = "ALA"),
              size = 0.5) +
  geom_violin(aes(x = 4,
                  y = C_20.5n.3,
                  fill = "EPA"),
              size = 0.5) +
  geom_violin(aes(x = 5,
                  y = C_22.6n.3,
                  fill = "DHA"),
              size = 0.25) +
  labs(x = "",
       y = "Proportion FA") +
  scale_fill_manual(values = colours,
                    name = "") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

dev.off()


