# Load libraries, data, etc. ----

library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(dplyr)
library(MuMIn)
library(tidyr)
library(boral)
library(R2jags)
library(coda)
library(lattice)
library(ggridges)

# Load data ----

FAs <- read.csv("FAs_concentration.csv", header = TRUE)
str(FAs)
FAs$ID <- as.factor(FAs$ID)
FAs$corral <- as.factor(FAs$corral)
FAs$sample.type <- as.factor(FAs$sample.type)
FAs$date <- as.Date(FAs$date,
                    format = "%b. %d, %Y")

treatments <- data.frame(corral = 
                           as.factor(c("B", "C", "D", "E", "F", "G", "H", "I")),
                         MPconcentration = as.numeric(c(0, 414, 29240, 100, 6, 
                                                        7071, 0, 1710)))

FAs <- left_join(FAs, 
                 treatments,
                 by = "corral")

# Assign corral A MP concentration of 24 for zoop_FA

FAs$MPconcentration[is.na(FAs$MPconcentration)] <- 24

# Separate out by sample type ----

perch_FA <- subset(FAs, sample.type == "perch")
zoop_FA <- subset(FAs, sample.type == "zooplankton")
food_FA <- subset(FAs, sample.type == "fish food")

# Add in perch biometrics data ----

perch_biometrics <- read.csv("perch_biometrics.csv", header = TRUE)

str(perch_biometrics)

perch_biometrics$corral <- as.factor(perch_biometrics$corral)
perch_biometrics$ID <- as.factor(perch_biometrics$ID)
perch_biometrics$sex <- as.factor(perch_biometrics$sex)

# Combine with perch FA data

perch_FA2 <- left_join(perch_FA, perch_biometrics, 
                       by = c("ID", "corral"))

perch_FA2$corral <- as.character(perch_FA2$corral)
perch_FA2$corral <- as.factor(perch_FA2$corral)

# Center and scale predictors

perch_FA2$scaled.date <- 
  as.numeric(scale(as.numeric(perch_FA2$date), 
                   center = min(as.numeric(perch_FA2$date)),
                   scale = 9))

perch_FA2$scaled.body.weight <-
  as.numeric(scale(perch_FA2$body.weight,
                   center = TRUE))

perch_FA2$scaled.MPconcentration <-
  as.numeric(scale(perch_FA2$MPconcentration,
                   center = TRUE))

Y_perch_FA <- 
  perch_FA2[,c(7:16,18:25,27:34,36:44)]  # pull out FA concentrations

X_perch_FA <- 
  perch_FA2[,c(60:62)]  # pull out predictors

ranef_ids_perch_FA <-
  as.integer(perch_FA2[,2])
  
offset_perch_FA <- 
  cbind(replicate(ncol(Y_perch_FA), perch_FA2[,50]))

boral.FA.mod <-
  boral(y = Y_perch_FA,
        X = X_perch_FA,,
        family = "normal",
        ranef.ids = ranef_ids_perch_FA,
        offset = offset_perch_FA,
        save.model = TRUE,
        mcmc.control = list(n.burnin = 1000,
                            n.iteration = 5000,
                            n.thin = 4,
                            seed = 31254))

## Specify model ####

perchFAmod <-
  function() {
    ## Data Level ##
    for (i in 1:n) {
      for (j in 1:p) {
        eta[i, j] <-
          ranef.coefs.ID1[j, ranef.ids[i]] +
          beta.date[j] * date[i] +
          beta.weight[j] * weight[i] +
          beta.conc[j] * conc[i] +
          offset[i, j]
        
        y[i, j] ~ dnorm(lv.coefs[j] + eta[i, j],
                        pow(lv.sd[j],-2))
        
      }
    }
    
    ## Process level and priors ##
    for (j in 1:p) {
      lv.coefs[j] ~ dnorm(0, 0.1)
      for (i in 1:n.ranefID) {
        ranef.coefs.ID1[j, i] ~ dnorm(0, pow(ranef.sigma.ID1[j],-2))
      }
      ranef.sigma.ID1[j] ~ dunif(0, 30)
      beta.date[j] ~ dnorm(0, 0.1)
      beta.weight[j] ~ dnorm(0, 0.1)
      beta.conc[j] ~ dnorm(0, 0.1)
      lv.sd[j] ~ dunif(0, 30) ## Dispersion parameters
    } ## Separate response intercepts
  }

## Generate initial values for MCMC ####

perchFAmodinit <- function()
{
  list(
    "lv.coefs" = rnorm(35),
    "ranef.sigma.ID1" = runif(35, min = 0, max = 30),
    "beta.date" = rnorm(35),
    "beta.weight" = rnorm(35),
    "beta.conc" = rnorm(35),
    "lv.sd" = runif(35, min = 0, max = 30)
  )
}

## Keep track of parameters ####

perchFAmodparam <- c("beta.date",
                     "beta.weight",
                     "beta.conc")

## Specify data for model ####

perchFAmoddata <-
  list(
    n = nrow(perch_FA2),
    p = ncol(Y_perch_FA),
    date = X_perch_FA$scaled.date,
    weight = X_perch_FA$scaled.body.weight,
    conc = X_perch_FA$scaled.MPconcentration,
    offset = offset_perch_FA,
    y = Y_perch_FA,
    n.ranefID = max(ranef_ids_perch_FA),
    ranef.ids = ranef_ids_perch_FA
  )

## Run the model ####

perchFAmodrun1 <-
  jags.parallel(
    data = perchFAmoddata,
    inits = perchFAmodinit,
    parameters.to.save = perchFAmodparam,
    n.chains = 3,
    n.cluster = 16,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 4,
    jags.seed = 46956,
    model = perchFAmod
    )

perchFAmodrun1
perchFAmodrun1mcmc <- as.mcmc(perchFAmodrun1)
xyplot(perchFAmodrun1mcmc)  # trace plots

perchFAmodrun1$BUGSoutput$sims.list