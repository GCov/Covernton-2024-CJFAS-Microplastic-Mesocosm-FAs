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
library(reshape2)

extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance", ]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}  # handy function for extracting posterior estimates of parameters

summarize.post <- function(x){
  x %>% 
    group_by(variable) %>% 
    summarize(mean = mean(value),
              lower95 = quantile(value, 0.025),
              upper95 = quantile(value, 0.979))
}

FA.names <- c("12:0", "13:0", "14:0", "15:0", "16:0", "17:0", "18:0", "20:0",
              "22:0", "24:0", "12:1", "14:1", "16:1n-7", "16:1n-9", "18:1n-7",
              "18:1n-9", "20:1n-9", "22:1n-9", "18:2n-6", "18:3n-6", "20:2n-6",
              "20:3n-6", "20:4n-6", "22:2n-6", "22:4n-6", "22:5n-6", "18:3n-3",
              "18:4n-3", "20:3n-3", "20:4n-3", "20:5n-3", "22:5n-3", "22:6n-3",
              "24:5n-3", "24:6n-3")

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

# Perch Model ----

## Prepare Data ----

Y_perch_FA <- 
  FAs_percent[,c(7:16,18:25,27:34,36:44)]  # pull out FA concentrations

X_perch_FA <- 
  perch_FA2[,c(61:62)]  # pull out predictors

ranef_ids_perch_FA <-
  as.integer(perch_FA2[,2])
  
offset_perch_FA <- 
  cbind(replicate(ncol(Y_perch_FA), perch_FA2[,50]))

## Specify model ####

perchFAmod <-
  function() {
    ## Data Level ##
    for (i in 1:n) {
      for (j in 1:p) {
        eta[i, j] <-
          ranef.coefs.ID1[j, ranef.ids[i]] +
          beta.weight[j] * weight[i] +
          beta.conc[j] * conc[i]
        
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
    "beta.weight" = rnorm(35),
    "beta.conc" = rnorm(35),
    "lv.sd" = runif(35, min = 0, max = 30)
  )
}

## Keep track of parameters ####

perchFAmodparam <- c("beta.weight",
                     "beta.conc")

## Specify data for model ####

perchFAmoddata <-
  list(
    n = nrow(perch_FA2),
    p = ncol(Y_perch_FA),
    weight = X_perch_FA$scaled.body.weight,
    conc = X_perch_FA$scaled.MPconcentration,
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
    n.iter = 40000,
    n.burnin = 1000,
    n.thin = 39,
    jags.seed = 46956,
    model = perchFAmod
    )

perchFAmodrun1
perchFAmodrun1mcmc <- as.mcmc(perchFAmodrun1)

perchFAmodrun1_post <- extract.post(perchFAmodrun1)

perchFAmodrun1_post_summary <- summarize.post(perchFAmodrun1_post)
  
perchFAmodrun1_post_summary$FA <- rep(FA.names, times = 2)

perchFAmodrun1_post_summary$coefficient <- 
  c(rep("MP Concentration", times = 35),
    rep("Body Weight", times = 35))

png("Perch Multivariate GLMM Plot.png",
    width = 19,
    height= 20, 
    units = "cm",
    res = 600)

ggplot(perchFAmodrun1_post_summary) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed") +
  geom_point(aes(x = mean, 
                 y = coefficient),
             colour = "red",
             size = 1.5) +
  geom_linerange(aes(xmin = lower95,
                     xmax = upper95,
                     y = coefficient),
                 colour = "red",
                 size = 0.5) +
  labs(x = "Slope",
       y = "Variable") +
  facet_wrap(~FA,
             ncol = 5,
             scales = "free_x") +
  theme1 +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1))
  

dev.off()

