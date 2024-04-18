library(ggplot2)
library(dplyr)
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

# Perch Model ----

# Center and scale predictors

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
