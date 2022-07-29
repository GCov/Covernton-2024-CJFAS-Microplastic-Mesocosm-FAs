# Load libraries, data, etc. ----

library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(Hmsc)
library(beepr)
library(lattice)
library(R2jags)
library(vegan)

set.seed(632)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  )

extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance" &
                 long$variable != "r", ]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}  # handy function for extracting posterior estimates of parameters


## Load data ----

perch_diet <- read.csv("perch_diet.csv", header = TRUE)

perch_diet$MPconcentration <- as.numeric(perch_diet$MPconcentration)
perch_diet$corral <- as.factor(perch_diet$corral)

# Count plots ----

## Plot totals

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = total.animals)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = cladocera)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

ggplot(data = perch_diet,
       aes(x = MPconcentration,
           y = cyclopoida)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1

## Put data into long format

perch_diet_long <- 
  perch_diet %>%
  pivot_longer(names(perch_diet[12:17]),
               names_to = "taxa",
               values_to = "count")

str(perch_diet_long)

perch_diet_long$taxa <- as.factor(perch_diet_long$taxa)
perch_diet_long$ID <- as.factor(perch_diet_long$ID)
perch_diet_long$count <- as.numeric(perch_diet_long$count)

levels(perch_diet_long$taxa) <-
  c("Amphipods",
    "Chironomid Larvae",
    "Chironomid Pupae",
    "Cladocerans",
    "Cyclopoid Copepods",
    "Larval Odonates")

perch_diet_long$treatment <-
  as.factor(perch_diet_long$corral)

levels(perch_diet_long$treatment) <-
  c("0(B)", "414", "29,240", "100", "6", "7,071", "0(H)", "1,710")

## Plot by taxa ----

png("Perch Diet Plot by Taxa.png",
    width = 19,
    height= 10, 
    units = "cm",
    res = 600)

ggplot(perch_diet_long) +
  geom_col(aes(x = ID,
               y = count,
               fill = reorder(taxa, 1/(count+1), mean)),
           colour = "black",
           size = 0.25) +
  labs(x = expression(paste("MP exposure concentration (particles"~L^-1*")")),
       y = "Number of Individuals (log scale)") +
  scale_fill_brewer(type = "qual",
                    name = "Taxa",
                    palette = "Set1") +
  scale_y_continuous(expand = c(0.005,0.005),
                     breaks = c(0, 1, 10, 100, 1000)) +
  coord_trans(y = "log1p") +
  facet_grid(.~reorder(treatment, MPconcentration, mean),
             scales = "free_x",
             switch = "x",
             space = "free_x") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

dev.off()


## Plot by relative abundance ----

perch_relabund <- perch_diet
perch_relabund[,c(12:17)] <- 
  decostand(perch_relabund[,c(12:17)], 
            method = "total")

perch_relabund_long <- 
  perch_relabund %>%
  pivot_longer(names(perch_relabund[12:17]),
               names_to = "taxa",
               values_to = "count")

perch_relabund_long$taxa <- as.factor(perch_relabund_long$taxa)
perch_relabund_long$ID <- as.factor(perch_relabund_long$ID)
perch_relabund_long$count <- as.numeric(perch_relabund_long$count)

levels(perch_relabund_long$taxa) <-
  c("Amphipods",
    "Chironomid Larvae",
    "Chironomid Pupae",
    "Cladocerans",
    "Cyclopoid Copepods",
    "Larval Odonates")

perch_relabund_long$treatment <-
  as.factor(perch_relabund_long$corral)

levels(perch_relabund_long$treatment) <-
  c("0(B)", "414", "29,240", "100", "6", "7,071", "0(H)", "1,710")

png("Perch Diet Plot by Taxa Relative Abundance.png",
    width = 33,
    height= 13, 
    units = "cm",
    res = 600)

ggplot(perch_relabund_long) +
  geom_col(aes(x = ID,
               y = count,
               fill = reorder(taxa, 1/(count+1), mean)),
           colour = "black",
           size = 0.25) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = expression(paste("Dose (MPs"~L^-1*")")),
       y = "Proportion of Individuals") +
  scale_fill_brewer(type = "qual",
                    name = "Taxa",
                    palette = "Set1") +
  facet_grid(.~reorder(treatment, dose, mean),
             scales = "free_x",
             switch = "x",
             space = "free_x") +
  theme1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

dev.off()






# GLMM for total individuals ----

# Standarize predictors
# Poisson
mod1 <- glmmTMB(total.animals ~ 
                  scale(MPconcentration, center = TRUE) + 
                  scale(body.length, center = TRUE) +  
                          (1 | corral),
                family = poisson(link = "log"),
                data = perch_diet)  
summary(mod1)
res1 <- simulateResiduals(mod1)
plot(res1)

# NB with linear variance
mod2 <- glmmTMB(total.animals ~ 
                  scale(MPconcentration, center = TRUE) + 
                  scale(body.length, center = TRUE) +  
                  (1 | corral),
                family = nbinom1(link = "log"),
                data = perch_diet)
summary(mod2)
res2 <- simulateResiduals(mod2)
plot(res2)

# NB with quadratic variance
mod3 <- glmmTMB(total.animals ~ 
                  scale(MPconcentration, center = TRUE) + 
                  scale(body.length, center = TRUE) +  
                  (1 | corral),
                family = nbinom2(link = "log"),
                data = perch_diet)
summary(mod3)
res3 <- simulateResiduals(mod3)
plot(res3)

AICc(mod2, mod3)  # mod3 is the better fit

# NB quadratic

mod4 <- glmmTMB(total.animals ~ 
                  scale(MPconcentration, center = TRUE) + 
                  I(scale(MPconcentration, center = TRUE)^2) +
                  scale(body.length, center = TRUE) +  
                  (1 | corral),
                family = nbinom2(link = "log"),
                data = perch_diet)

summary(mod4)
res4 <- simulateResiduals(mod4)
plot(res4)

AICc(mod3, mod4)  # quadratic model not a better fit

## Plot model predictions ----

simdata <- data.frame(MPconcentration = seq(0, max(perch_diet$MPconcentration), 
                                            length.out = 1000),
                      body.length = rep(mean(perch_diet$body.length), 1000),
                      corral = rep(NA, 1000))

simdata$mean <- predict(mod3, simdata)

simdata$se <- as.numeric(predict(mod3,simdata, se.fit = TRUE)$se.fit)

simdata$upper <- with(simdata, mean + 1.96*se)
simdata$lower <- with(simdata, mean - 1.96*se)

png("Perch Diet Totals.png",
    width = 9,
    height= 7, 
    units = "cm",
    res = 600)

ggplot() +
  geom_ribbon(data = simdata,
              aes(x = MPconcentration,
                  ymin = exp(lower),
                  ymax = exp(upper)),
              fill = "lime green",
              alpha = 0.3) +
  geom_line(data = simdata,
            aes(x = MPconcentration,
                y = exp(mean)),
            size = 0.5,
            linetype = "dashed") +
  geom_point(data = perch_diet,
             aes(x = MPconcentration,
                 y = total.animals)) +
  scale_y_continuous(expand = c(0.009,0.005),
                     breaks = c(0, 1, 10, 100, 1000)) +
  coord_trans(y = "log1p") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000),
                     expand = c(0.01, 0.01)) +
  labs(x = expression(paste("MP Exposure Concentration (particles"~L^-1*")")),
       y = "Total number of Individuals (log scale)") +
  theme1 +
  theme(plot.margin = margin(0.1, 0.7, 0.1, 0.1, unit = "cm"))

dev.off()

# Check that total individuals is a good measure of gut fullness ----

gutmod1 <- lm(log(gi.weight) ~ log(total.animals + 1), data = perch_diet)

plot(resid(gutmod1) ~ fitted(gutmod1), data = perch_diet)
abline(0,0)

gutmodresid <- simulateResiduals(gutmod1)
plot(gutmodresid)  # looks good

summary(gutmod1)  
# total indv. is a pretty good predictor of GI weight (R2 = 0.38)

gutpred <- predict(gutmod1, se.fit = TRUE)

gutpred$lower <- with(gutpred, fit - se.fit)
gutpred$upper <- with(gutpred, fit + se.fit)


## Plot model predictions ----

png("Fullness Plot.png",
    width = 17,
    height= 13, 
    units = "cm",
    res = 600)

ggplot(data = perch_diet) +
  geom_point(aes(x = total.animals,
                 y = gi.weight)) +
  geom_ribbon(aes(x = total.animals,
                  ymin = exp(gutpred$lower),
                  ymax = exp(gutpred$upper)),
              fill = "aquamarine",
              alpha = 0.3) +
  geom_line(aes(x = total.animals,
                y = exp(gutpred$fit)),
            size = 0.5,
            linetype = "dashed") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  labs(x = "Total Number of Individuals (log scale)",
       y = "Gastrointestinal Tract Weight (g)") +
  theme1

dev.off()



# GLMM by taxa ----

## Poisson model ----

### Fit model ----

taxmod1 <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- 
      alpha_taxa[taxa[i]] +
      beta_dose[taxa[i]] * dose[i] +
      beta_length * length[i] +
      gamma_ind[ind[i]]
    
    gamma_ind[i] ~ dnorm(gamma_corral[corral[i]], tau_ind)
    
    # Fitted values
    fitted[i] ~ dpois(lambda[i])
  }
  
  # Priors
  
  for(j in 1:ntaxa) {
    alpha_taxa[j] ~ dnorm(0, 1)
    beta_dose[j] ~ dnorm(0, 1)
  }
  
  beta_length ~ dnorm(0, 1)
  
  for(k in 1:ncorral) {
    gamma_corral[k] ~ dnorm(0, tau_corral)
  }
  
  tau_ind <- inverse(pow(sigma_ind, 2))
  sigma_ind ~ dexp(1)
  
  tau_corral <- inverse(pow(sigma_corral, 2))
  sigma_corral ~ dexp(1)
}

### Generate initial values for MCMC ----

taxmod1init <- function()
{
  list(
    "alpha_taxa" = rnorm(6),
    "beta_dose" = rnorm(6),
    "beta_length" = rnorm(1),
    "sigma_ind" = rexp(1),
    "sigma_corral" = rexp(1)
  )
}

### Keep track of parameters ----

taxmod1param <- c("alpha_taxa",
                  "beta_dose",
                  "beta_length",
                  "sigma_ind",
                  "sigma_corral")

### Specify data ----

taxmod1data <-
  list(
    y = perch_diet_long$count,
    N = nrow(perch_diet_long),
    taxa = as.integer(perch_diet_long$taxa),
    dose = as.numeric(scale(perch_diet_long$MPconcentration, center = TRUE)),
    length = as.numeric(scale(perch_diet_long$body.length, center = TRUE)),
    ind = as.integer(perch_diet_long$ID),
    corral = as.integer(perch_diet_long$corral),
    ntaxa = length(unique(perch_diet_long$taxa)),
    ncorral = length(unique(perch_diet_long$corral))
  )

### Run the model ----
taxmod1run1 <- jags.parallel(
  data = taxmod1data,
  inits = taxmod1init,
  parameters.to.save = taxmod1param,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 32165,
  model = taxmod1
)

taxmod1run1
taxmod1run1mcmc <- as.mcmc(taxmod1run1)
xyplot(taxmod1run1mcmc, 
       layout = c(6, ceiling(nvar(taxmod1run1mcmc)/6)))  # trace plots

### Diagnostics ----

taxmod1param2 <- c("fitted", "lambda")

taxmod1run2 <- jags.parallel(
  data = taxmod1data,
  inits = taxmod1init,
  parameters.to.save = taxmod1param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 2000,
  n.burnin = 500,
  n.thin = 1,
  jags.seed = 32165,
  model = taxmod1
)

taxmod1.response <- t(taxmod1run2$BUGSoutput$sims.list$fitted)
taxmod1.observed <- perch_diet_long$count
taxmod1.fitted <- apply(t(taxmod1run2$BUGSoutput$sims.list$lambda),
                        1,
                        median)

check.taxmod1 <- createDHARMa(
  simulatedResponse = taxmod1.response,
  observedResponse = taxmod1.observed,
  fittedPredictedResponse = taxmod1.fitted,
  integerResponse = T
)

plot(check.taxmod1)  # looks pretty bad.. zero-inflated Poisson, NB?
testZeroInflation(check.taxmod1)
testDispersion(check.taxmod1)


## Negative binomial model ----

### Fit model ----

taxmod2 <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r/(r+lambda[i])
    log(lambda[i]) <- 
      alpha_taxa[taxa[i]] +
      beta_dose[taxa[i]] * dose[i] +
      beta_length * length[i] +
      gamma_ind[ind[i]]
    
    gamma_ind[i] ~ dnorm(gamma_corral[corral[i]], tau_ind)
    
    ## Fitted values
    fitted[i] ~ dpois(lambda[i])
  }
  
  # Priors
  
  r ~ dunif(0,50)
  
  for(j in 1:ntaxa) {
    alpha_taxa[j] ~ dnorm(0, 1)
    beta_dose[j] ~ dnorm(0, 1)
  }
  
  beta_length ~ dnorm(0, 1)
  
  for(k in 1:ncorral) {
    gamma_corral[k] ~ dnorm(0, tau_corral)
  }
  
  tau_ind <- inverse(pow(sigma_ind, 2))
  sigma_ind ~ dexp(1)
  
  tau_corral <- inverse(pow(sigma_corral, 2))
  sigma_corral ~ dexp(1)
}

### Generate initial values for MCMC ----

taxmod2init <- function()
{
  list(
    "r" = runif(1, 0, 50),
    "alpha_taxa" = rnorm(6),
    "beta_dose" = rnorm(6),
    "beta_length" = rnorm(1),
    "sigma_ind" = rexp(1),
    "sigma_corral" = rexp(1)
  )
}

#### Keep track of parameters ----

taxmod2param <- c("alpha_taxa",
                  "beta_dose",
                  "beta_length",
                  "sigma_ind",
                  "sigma_corral")

#### Specify data ----

taxmod2data <-
  list(
    y = perch_diet_long$count,
    N = nrow(perch_diet_long),
    taxa = as.integer(perch_diet_long$taxa),
    dose = as.numeric(scale(perch_diet_long$MPconcentration, center = TRUE)),
    length = as.numeric(scale(perch_diet_long$body.length, center = TRUE)),
    ind = as.integer(perch_diet_long$ID),
    corral = as.integer(perch_diet_long$corral),
    ntaxa = length(unique(perch_diet_long$taxa)),
    ncorral = length(unique(perch_diet_long$corral))
  )

#### Run the model ----
taxmod2run1 <- jags.parallel(
  data = taxmod2data,
  inits = taxmod2init,
  parameters.to.save = taxmod2param,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 4,
  jags.seed = 32165,
  model = taxmod2
)

taxmod2run1
taxmod2run1mcmc <- as.mcmc(taxmod2run1)
xyplot(taxmod2run1mcmc, 
       layout = c(6, ceiling(nvar(taxmod2run1mcmc)/6)))  # trace plots

#### Diagnostics ----

taxmod2param2 <- c("fitted", "lambda")

taxmod2run2 <- jags.parallel(
  data = taxmod2data,
  inits = taxmod2init,
  parameters.to.save = taxmod2param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 5000,
  n.burnin = 500,
  n.thin = 4,
  jags.seed = 32165,
  model = taxmod2
)

taxmod2.response <- t(taxmod2run2$BUGSoutput$sims.list$fitted)
taxmod2.observed <- perch_diet_long$count
taxmod2.fitted <- apply(t(taxmod2run2$BUGSoutput$sims.list$lambda),
                        1,
                        median)

check.taxmod2 <- createDHARMa(
  simulatedResponse = taxmod2.response,
  observedResponse = taxmod2.observed,
  fittedPredictedResponse = taxmod2.fitted,
  integerResponse = T
)

plot(check.taxmod2)
testZeroInflation(check.taxmod2)  ## zero-inflated


## Zero-inflated Poisson ----

### Fit model ----

taxmod3 <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(mu[i])
    mu[i] <- lambda[i] * (1 - z[i]) + 0.0000000001
    z[i] ~ dbern(p[i])
    logit(p[i]) <- alpha_p
    log(lambda[i]) <- 
      alpha_taxa[taxa[i]] +
      beta_dose[taxa[i]] * dose[i] +
      beta_length * length[i] +
      gamma_ind[ind[i]]
    
    gamma_ind[i] ~ dnorm(gamma_corral[corral[i]], tau_ind)
    
    ## Fitted values
    fitted[i] ~ dpois(mu[i])
  }
  
  # Priors
  
  alpha_p ~ dnorm(0, 0.2)
  
  for(j in 1:ntaxa) {
    alpha_taxa[j] ~ dnorm(0, 1)
    beta_dose[j] ~ dnorm(0, 1)
  }
  
  beta_length ~ dnorm(0, 1)
  
  for(k in 1:ncorral) {
    gamma_corral[k] ~ dnorm(0, tau_corral)
  }
  
  tau_ind <- inverse(pow(sigma_ind, 2))
  sigma_ind ~ dexp(1)
  
  tau_corral <- inverse(pow(sigma_corral, 2))
  sigma_corral ~ dexp(1)
}

### Generate initial values for MCMC ----

taxmod3init <- function()
{
  list(
    "alpha_p" = rnorm(1, 0, 5),
    "alpha_taxa" = rnorm(6),
    "beta_dose" = rnorm(6),
    "beta_length" = rnorm(1),
    "sigma_ind" = rexp(1),
    "sigma_corral" = rexp(1)
  )
}

### Keep track of parameters ----

taxmod3param <- c("alpha_p", "
                  alpha_taxa",
                  "beta_dose",
                  "beta_length",
                  "sigma_ind",
                  "sigma_corral")

### Specify data ----

taxmod3data <-
  list(
    y = perch_diet_long$count,
    N = nrow(perch_diet_long),
    taxa = as.integer(perch_diet_long$taxa),
    dose = as.numeric(scale(perch_diet_long$MPconcentration, center = TRUE)),
    length = as.numeric(scale(perch_diet_long$body.length, center = TRUE)),
    ind = as.integer(perch_diet_long$ID),
    corral = as.integer(perch_diet_long$corral),
    ntaxa = length(unique(perch_diet_long$taxa)),
    ncorral = length(unique(perch_diet_long$corral))
  )

### Run the model ----
taxmod3run1 <- jags.parallel(
  data = taxmod3data,
  inits = taxmod3init,
  parameters.to.save = taxmod3param,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 100000,
  n.burnin = 5000,
  n.thin = 95,
  jags.seed = 32165,
  model = taxmod3
)

taxmod3run1
taxmod3run1mcmc <- as.mcmc(taxmod3run1)
xyplot(taxmod3run1mcmc, 
       layout = c(6, ceiling(nvar(taxmod3run1mcmc)/6)))  # trace plots

### Diagnostics ----

taxmod3param2 <- c("fitted", "mu")

taxmod3run2 <- jags.parallel(
  data = taxmod3data,
  inits = taxmod3init,
  parameters.to.save = taxmod3param2,
  n.chains = 3,
  n.cluster = 16,
  n.iter = 100000,
  n.burnin = 5000,
  n.thin = 95,
  jags.seed = 32165,
  model = taxmod3
)

taxmod3.response <- t(taxmod3run2$BUGSoutput$sims.list$fitted)
taxmod3.observed <- perch_diet_long$count
taxmod3.fitted <- apply(t(taxmod3run2$BUGSoutput$sims.list$mu),
                        1,
                        median)

check.taxmod3 <- createDHARMa(
  simulatedResponse = taxmod3.response,
  observedResponse = taxmod3.observed,
  fittedPredictedResponse = taxmod3.fitted,
  integerResponse = T
)

plot(check.taxmod3)  ## also sucks




# Try HMSC model ----

## Specify data ----

XData <- perch_diet[,c(3,5)]  # pull out predictors

design <- data.frame(corral = factor(perch_diet[,2]))
rlevels <- HmscRandomLevel(units = design$corral)

Y <- perch_diet[,c(12:17)]  # pull out species counts

## Specify model structure ----

model1 <- Hmsc(Y = Y,  # response data
               XData = XData,  # covariates
               XFormula = ~dose + body.length,  # model formula
               XScale = TRUE,  # scale covariates for fixed effects,
               studyDesign = design,
               ranLevels = list(corral = rlevels),
               distr = "poisson")

## Run MCMC chains ----

run1 <- sampleMcmc(model1,
                   thin = 20,
                   samples = 20000,
                   transient = 500,
                   nChains = 3,
                   nParallel = 3,
                   verbose = 1000)
beep(8)

## Check convergence ----

post1 <- convertToCodaObject(run1)

effectiveSize(post1$Beta)
gelman.diag(post1$Beta, 
            transform = TRUE,
            multivariate = FALSE)$psrf

plot(post1$Beta)


## Assess explanatory power ----

pred1 <- computePredictedValues(run1)
evaluateModelFit(hM = run1, predY = pred1)


## Cross validation ----

partition1 <- createPartition(run1, nfolds = 2)
pred1.1 <- computePredictedValues(run1, partition = partition1)

evaluateModelFit(hM = run1, predY = pred1.1)  # predictive power sucks

## Look at slope estimates ----
 
postBeta <- getPostEstimate(run1, parName = "Beta")
plotBeta(run1, post = postBeta, param = "Support", supportLevel = 0.95)


# Look at effect of body size ----

ggplot(perch_diet_long) +
  geom_point(aes(x = body.length,
                 y = count,
                 fill = reorder(taxa, 1/(count+1), mean)),
             colour = "black",
             shape = 21,
             size = 4) +
  labs(x = "Total Length (cm)",
       y = "Number of Individuals (log scale)") +
  scale_fill_brewer(type = "qual",
                    name = "Taxa",
                    palette = "Set1") +
  theme1 +
  theme(legend.position = "bottom")



# nMDS ----

# Remove rows with all zeros
sum <- numeric()
for(i in 1:nrow(Y)){
  sum[i] <- sum(Y[i,])
}

Y2 <- Y[sum>0,]
XData2 <- XData[sum>0,]

# Convert to relative abundance
rel_abund <- decostand(Y2, method = "total")

# Calculate distance matrix
diss <- as.matrix(vegdist(rel_abund, method = "bray", na.rm = TRUE), 
                  labels = TRUE)

nMDS1 <- metaMDS(diss,
                 distance = "bray",
                 k = 3,
                 maxis = 999,
                 trymax = 250,
                 wascores = TRUE)

# Shepards test/goodness of fit
goodness(nMDS1)
stressplot(nMDS1)

data.scores <- as.data.frame(scores(nMDS1))
data.scores$dose <- XData2$dose
data.scores$body.length <- XData2$body.length

## Generate hulls ----

data.scores$dose <- as.factor(data.scores$dose)

hulls <- data.frame()

for(i in 1:length(unique(data.scores$dose))) {
  hull <- 
    data.scores[data.scores$dose ==
                  unique(data.scores$dose)[i],
                ][chull(data.scores[data.scores$dose ==
                                      unique(data.scores$dose)[i],
                                    c(1:2)]), ]
  hulls <- rbind(hulls, hull)
}

## Plot ----

png("Perch nMDS.png",
    width = 20,
    height= 18, 
    units = "cm",
    res = 600)

ggplot() +
  geom_polygon(data = hulls,
               aes(x = NMDS1,
                   y = NMDS2,
                   fill = dose,
                   colour = dose),
               alpha = 0.3,
               size = 0.5) +
  geom_label(data = data.scores,
             aes(x = NMDS1,
                 y = NMDS2,
                 label = dose,
                 fill = dose,
                 size = body.length),
             alpha = 0.3) +
  scale_size(name = "Total Length (cm)") +
  scale_fill_brewer(type = "seq",
                    palette = "YlOrRd",
                    name = expression(paste("Dose (MPs"~L^-1*")"))) +
  scale_colour_brewer(type = "seq",
                      palette = "YlOrRd",
                      name = expression(paste("Dose (MPs"~L^-1*")"))) +
  theme1

dev.off()


