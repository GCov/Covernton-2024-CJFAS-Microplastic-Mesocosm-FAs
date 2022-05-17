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
library(Hmsc)

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

FA.names <- c("12:0", 
              "24:0", "12:1", "14:1", "16:1n-7", "16:1n-9", "18:1n-7",
              "18:1n-9", "20:1n-9", "22:1n-9", "18:2n-6", 
              "13:0", 
              "18:3n-6", "20:2n-6",
              "20:3n-6", "20:4n-6", "22:2n-6", "22:4n-6", "22:5n-6", "18:3n-3",
              "18:4n-3", "20:3n-3", 
              "14:0", 
              "20:4n-3", "20:5n-3", "22:5n-3", "22:6n-3",
              "24:5n-3", "24:6n-3",
              "15:0", "16:0", "17:0", "18:0", "20:0", "22:0")

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
  perch_FA2[,c(7:16,18:25,27:34,36:44)]  # pull out FA concentrations

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
          intercept[j] +
          beta.weight[j] * weight[i] +
          beta.conc[j] * conc[i] +
          ranef.coefs.ID1[j, ranef.ids[i]]
        y[i, j] ~ dnorm(lv.coefs[j] + eta[i, j],
                        pow(lv.sd[j],-2))
        
      }
    }
    
    ## Process level and priors ##
    for (j in 1:p) {
      intercept[j] ~ dnorm(0, 1)
      lv.coefs[j] ~ dnorm(0, 1)
      for (i in 1:n.ranefID) {
        ranef.coefs.ID1[j, i] ~ dnorm(0, pow(ranef.sigma.ID1[j],-2))
      }
      ranef.sigma.ID1[j] ~ dunif(0, 1)
      beta.weight[j] ~ dnorm(0, 1)
      beta.conc[j] ~ dnorm(0, 1)
      lv.sd[j] ~ dunif(0, 1) ## Dispersion parameters
    } ## Separate response intercepts
  }

## Generate initial values for MCMC ####

perchFAmodinit <- function()
{
  list(
    "intercept" = rnorm(35),
    "lv.coefs" = rnorm(35),
    "ranef.sigma.ID1" = runif(35),
    "beta.weight" = rnorm(35),
    "beta.conc" = rnorm(35),
    "lv.sd" = runif(35)
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
    y = log(Y_perch_FA + 1),
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

xyplot(perchFAmodrun1mcmc, 
       layout = c(6, ceiling(nvar(perchFAmodrun1mcmc)/6)))

perchFAmodrun1_post <- extract.post(perchFAmodrun1)

perchFAmodrun1_post_summary <- summarize.post(perchFAmodrun1_post)
  
perchFAmodrun1_post_summary$FA <- rep(FA.names, 
                                      times = length(perchFAmodparam))

## Plot parameter posterior estimates ----

perchFAmodrun1_post_summary$coefficient <- 
  c(rep("MP Concentration", times = 35),
    rep("Body Weight", times = 35))

png("Perch Multivariate GLMM Parameter Posteriors Plot.png",
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
  labs(x = "Parameter Estimate",
       y = "Variable") +
  facet_wrap(~FA,
             ncol = 5,
             scales = "free_x") +
  theme1
  

dev.off()

## Predict from model ----

perch_FA_conc_new <- perch_FA2[,c(1:6, 51:62)]

perch_FA_conc_new$scaled.body.weight <- 
  mean(perch_FA_conc_new$scaled.body.weight)

# Pull out all necessary terms

perchFAmodparam2 <- c("beta.weight",
                     "beta.conc",
                     "intercept")

perchFAmodrun2 <-  # rerun model
  jags.parallel(
    data = perchFAmoddata,
    inits = perchFAmodinit,
    parameters.to.save = perchFAmodparam2,
    n.chains = 3,
    n.cluster = 16,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 4,
    jags.seed = 46956,
    model = perchFAmod
  )

perch_FA_conc_predictions <- list(mean = data.frame(),
                                  upper = data.frame(),
                                  lower = data.frame())

# Use posteriors to predict

for(i in 1:nrow(perch_FA_conc_new)) {
  for(j in 1:ncol(Y_perch_FA)){
  eta <- 
    perchFAmodrun2$BUGSoutput$sims.list$intercept[,j] +
    perchFAmodrun2$BUGSoutput$sims.list$beta.weight[,j] * 
    perch_FA_conc_new$scaled.body.weight[i] +
    perchFAmodrun2$BUGSoutput$sims.list$beta.conc[,j] * 
    perch_FA_conc_new$scaled.MPconcentration[i]
  perch_FA_conc_predictions$mean[i, j] <- exp(mean(eta)) - 1
  perch_FA_conc_predictions$upper[i, j] <- exp(quantile(eta, 0.975)) - 1
  perch_FA_conc_predictions$lower[i, j] <- exp(quantile(eta, 0.025)) - 1
}
}

# rename columns

colnames(perch_FA_conc_predictions$mean) <- colnames(Y_perch_FA)
colnames(perch_FA_conc_predictions$upper) <- colnames(Y_perch_FA)
colnames(perch_FA_conc_predictions$lower) <- colnames(Y_perch_FA)

# bind everything together

perch_FA_conc_new <- cbind(perch_FA_conc_new, perch_FA_conc_predictions$mean)
  
perch_FA_conc_new_long <- 
  perch_FA_conc_new %>%
  pivot_longer(colnames(Y_perch_FA),
               names_to = "FA",
               values_to = "Mean")

perch_FA_conc_multi_upper <-
  perch_FA_conc_predictions$upper %>% 
  pivot_longer(colnames(Y_perch_FA),
               names_to = "FA",
               values_to = "Upper")

perch_FA_conc_multi_lower <-
  perch_FA_conc_predictions$lower %>% 
  pivot_longer(colnames(Y_perch_FA),
               names_to = "FA",
               values_to = "Lower")

perch_FA_conc_new_long <- 
  cbind(perch_FA_conc_new_long, perch_FA_conc_multi_upper[,2])

perch_FA_conc_new_long <- 
  cbind(perch_FA_conc_new_long, perch_FA_conc_multi_lower[,2])

# Put original data into long form

perch_FA2_long <- 
  perch_FA2 %>%
  pivot_longer(colnames(Y_perch_FA),
               names_to = "FA",
               values_to = "Value")

## Plot predictions ----

png("Perch Multivariate GLMM Predictions Plot.png",
    width = 19,
    height= 19, 
    units = "cm",
    res = 600)

ggplot() +
  geom_ribbon(data = perch_FA_conc_new_long,
              aes(x = MPconcentration,
                  ymin = Lower,
                  ymax = Upper),
              fill = "red",
              alpha = 0.3) +
  geom_line(data = perch_FA_conc_new_long,
            aes(x = MPconcentration,
                y = Mean),
            colour = "red") +
  geom_point(data = perch_FA2_long,
             aes(x = MPconcentration,
                 y = Value)) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = expression(paste("Concentration (mg "~g^-1*")"))) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ FA, ncol = 5) +
  theme1

dev.off()


# Perch HMSC model----

## Specify data ----

perch_FA_conc_predictors <- 
  perch_FA2[,c(51,54)]  # pull out predictors

perch_FA_conc_RE <- data.frame(corral = factor(perch_FA2[,2]))
perch_FA_conc_rlevels <- HmscRandomLevel(units = perch_FA_conc_RE$corral)

perch_FA_conc_response <- 
  perch_FA2[,c(7:16,18:25,27:34,36:44)]  # pull out FA proportions

## Specify model structure ----

model1 <- 
  Hmsc(Y = perch_FA_conc_response,  # response data
       XData = perch_FA_conc_predictors,  # covariates
       XFormula = ~ body.weight + MPconcentration,  # model formula
       XScale = TRUE,  # scale covariates for fixed effects,
       studyDesign = perch_FA_conc_RE,
       ranLevels = list(corral = perch_FA_conc_rlevels),
       distr = "normal")

## Run MCMC chains ----

set.seed(6461)

perch_FA_conc_run1 <- sampleMcmc(model1,
                                 thin = 1,
                                 samples = 2000,
                                 transient = 100,
                                 nChains = 3,
                                 nParallel = 3,
                                 verbose = 1000)

## Check convergence ----

perch_FA_conc_post1 <- convertToCodaObject(perch_FA_prop_run1)

effectiveSize(perch_FA_conc_post1$Beta)
gelman.diag(perch_FA_conc_post1$Beta, 
            transform = TRUE,
            multivariate = FALSE)$psrf

plot(perch_FA_conc_post1perch_FA_prop_post1$Beta)


## Assess explanatory power ----

perch_FA_conc_pred1 <- computePredictedValues(perch_FA_conc_run1)
evaluateModelFit(hM = perch_FA_conc_run1, predY = perch_FA_conc_pred1)


## Cross validation ----

partition1 <- createPartition(perch_FA_conc_run1, nfolds = 2)
perch_FA_conc_pred1.1 <- 
  computePredictedValues(perch_FA_conc_run1, partition = partition1)

evaluateModelFit(hM = perch_FA_conc_run1, 
                 predY = perch_FA_conc_pred1.1)

## Look at slope estimates ----

perch_FA_conc_postBeta <- 
  getPostEstimate(perch_FA_conc_run1, parName = "Beta")
plotBeta(perch_FA_conc_run1, post = perch_FA_conc_postBeta, 
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



