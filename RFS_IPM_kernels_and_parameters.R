##################################################################################################################
# Modeling Reticulated Flatwoods Salamander populations: Assessing effects of future climate change scenarios
# on long-term population viability. This script has just the final kernels and parameters that are needed
# to run the model. See full code for calculations of certain values.
##################################################################################################################


#-----------------------------------------------------------------------------------------------------------------
# 1.0 Setup Workspace
#-----------------------------------------------------------------------------------------------------------------

# Load Packages
library(tidyverse)
library(lubridate)
library(patchwork)


# functions
std <- function(x){
  (x-mean(x)/sd(x))
}

# Outline
# 1.0 Setup Workspace - run every time!
# 2.0 Fecundity kernel
# 3.0 Growth kernel
# 4.0 Survival kernel
# 5.0 Combine kernels and parameters
# 6.0 Extra functions



#-----------------------------------------------------------------------------------------------------------------
# 2.0 Fecundity kernel - Eggs produced, survival until metamorphosis, size at metamorphosis
#-----------------------------------------------------------------------------------------------------------------

# 2.1 Linear relationship defining how many eggs are produced by female salamanders of reproductive size ---------

# This variability is approximately a 50% CI around the intercept
FUN.eggs <- function(z, m.par) {
  int <- runif(1, (m.par[["Eggs.int"]] - 0.228), (m.par[["Eggs.int"]] + 0.228))
  p <- ifelse(z < m.par[["SVL.fem.min"]],
              0, 
              exp(m.par[["Eggs.slope"]] * z + int))
  return(p)
}


# 2.2 Function that estimates how many eggs survive until metamorphosis ------------------------------------------

# Random draw between 0.1 and 0.5
FUN.prob.rec <- function(x) {
  runif(x, 0.01, 0.05)
}

# Fixed recruitment at 0.01 intervals from 0.01 to 0.05
FUN.prob.rec.01 <- 0.01
FUN.prob.rec.02 <- 0.02
FUN.prob.rec.03 <- 0.03
FUN.prob.rec.04 <- 0.04
FUN.prob.rec.05 <- 0.05

# 2.3A Maturity function ------------------------------

z = i.par$meshpoints
A = -70
B = 1.5
p.mat = 1/(1+exp(-(A + B*z)))
p.mat
plot(z, p.mat)

# 2.3B Do female salamander reproduce every year after they reach a sufficient size? ------------------------------

# Three levels at 100%, 75%, and 50%
FUN.prob.breed.100 <- function(z, m.par){
  Prob.breed = 1/(1+exp(-(rnorm(1, m.par[["A"]],2) + m.par[["B"]]*z)))
  return(as.numeric(Prob.breed))
}

FUN.prob.breed.75 <- function(z, m.par){
  Prob.breed = 1/(1+exp(-((rnorm(1, m.par[["A"]],2) + m.par[["B"]]*z)))) * 0.75
  return(as.numeric(Prob.breed))
}

FUN.prob.breed.50 <- function(z, m.par){
  Prob.breed = 1/(1+exp(-((rnorm(1, m.par[["A"]],2) + m.par[["B"]]*z)))) * 0.5
  return(as.numeric(Prob.breed))
}


# 2.4 Functions that estimate metamorph size based on hydroperiod ------------------------------------------------

# Functions with variability in mean recruitment size
FUN.metas.long <- function(z1, m.par) {
  mu <- runif(1, (m.par[["Est.longhydro.mean"]] - 2), (m.par[["Est.longhydro.mean"]] + 2)) # mean size of recruits
  sigma <- m.par[["Est.longhydro.sd"]] # sd about mean
  p.den.metas.long <- dnorm(z1, mean = mu, sd = sigma) # pdf for offspring size z1
  return(p.den.metas.long)
}

FUN.metas.short <- function(z1, m.par) {
  mu <- runif(1, (m.par[["Est.medhydro.mean"]] - 2), (m.par[["Est.medhydro.mean"]] + 2)) # mean size of recruits
  sigma <- m.par[["Est.medhydro.sd"]] # sd about mean
  p.den.metas.short <- dnorm(z1, mean = mu, sd = sigma) # pdf for offspring size z1
  return(p.den.metas.short)
}


# Functions with constant mean recruitment size
FUN.metas.long.det <- function(z1, m.par) {
  mu <- m.par[["Est.longhydro.mean"]] # mean size of recruits
  sigma <- m.par[["Est.longhydro.sd"]] # sd about mean
  p.den.metas.long <- dnorm(z1, mean = mu, sd = sigma) # pdf for offspring size z1
  return(p.den.metas.long)
}

FUN.metas.short.det <- function(z1, m.par) {
  mu <- m.par[["Est.medhydro.mean"]] # mean size of recruits
  sigma <- m.par[["Est.medhydro.sd"]] # sd about mean
  p.den.metas.short <- dnorm(z1, mean = mu, sd = sigma) # pdf for offspring size z1
  return(p.den.metas.short)
}

# 2.5 Combine above functions to create fecundity kernel ---------------------------------------------------------

# Functions that combine the above portions of the reproductive kernel. There are functions for both a 'short' and
# 'long' hydroperiod, which impacts the size at metamorphosis.

# For a given individual of size z the probability of producing an individual of size z1 equals: 
#       p(breeding) * n(eggs) * p(recruitment) * recruit size distribution
# The function calculated total number of expected metamorphs and multiplies that by the expected size distribution
# based on the hydroperiod of the breeding wetland.

# Define function for the potential effects of flooding mismatch
FUN.flood <- function(x) {
  runif(x, 0.2, 1.0)
}

# Fully stochastic functions
FUN.fec.z.short <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.short(z1, m.par))
}

FUN.fec.z.long <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.long(z1, m.par))
}

FUN.fec.z.short.flood <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.short(z1, m.par) * FUN.flood(1))
}

FUN.fec.z.long.flood <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.long(z1, m.par) * FUN.flood(1))
}

# Functions to look at differences in reproduction parameters (probability that females breed and probability of recruitment)

# 100% of females breed
FUN.fec.z.100.01.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.100.02.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.100.03.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.100.04.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.100.05.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.100.01.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.100.02.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.100.03.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.100.04.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.100.05.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.100(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.short(z1, m.par))
}

# 75% of females breed
FUN.fec.z.75.01.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.75.02.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.75.03.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.75.04.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.75.05.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.75.01.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.75.02.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.75.03.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.75.04.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.75.05.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.75(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.short(z1, m.par))
}

# 50% of females breed
FUN.fec.z.50.01.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.50.02.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.50.03.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.50.04.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.50.05.l <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.long(z1, m.par))
}

FUN.fec.z.50.01.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.01 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.50.02.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.02 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.50.03.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.03 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.50.04.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.04 * FUN.metas.short(z1, m.par))
}

FUN.fec.z.50.05.s <- function(z1, z, m.par){
  return(0.5 * FUN.prob.breed.50(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec.05 * FUN.metas.short(z1, m.par))
}


#*CAN ADD SURVIVAL TO THESE BUT THIS IS REALLY CAPTURED IN THE PROBABILITY TO RECRUIT FUNCTION**
#FUN.fec.z.short <- function(z1, z, m.par){
#  return((1/2) * FUN.prob.breed(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.short(z1, m.par) * FUN.surv(z, m.par))
#}

#FUN.fec.z.long <- function(z1, z, m.par){
#  return((1/2) * FUN.prob.breed(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.long(z1, m.par) * FUN.surv(z, m.par))
#}

# 2.6 Define parameters needed for the fecundity kernel ----------------------------------------------------------

Fec.par <- c(
  # Egg Production
  Eggs.slope = 0.050296,
  Eggs.int = 1.952206,
  # Is female large enough to reproduce?
  A = -70,
  B = 1.5,
  SVL.fem.min = 47.3,
  # Percent of eggs surviving to metamorphosis
  FUN.prob.rec = FUN.prob.rec,
  # Means of metamorph size distribution
  Est.longhydro.mean = 42.69894,
  Est.medhydro.mean = 35.09483,
  # Standard deviations of metamorph size distribution
  Est.longhydro.sd = 3.710759,
  Est.medhydro.sd = 2.737135)

names(Fec.par) <- c("Eggs.slope", "Eggs.int", "A", "B", "SVL.fem.min", "Prob.rec", "Est.longhydro.mean", "Est.medhydro.mean",
                    "Est.longhydro.sd", "Est.medhydro.sd")


#-----------------------------------------------------------------------------------------------------------------
# 3.0 Growth kernel - How do salamanders grow after metamorphosis?
#-----------------------------------------------------------------------------------------------------------------

# 3.1 Load posterior distribution of von Bertalanffy parameters (Brooks et al. 2020) -----------------------------
load("Data/RFS_Linf_post.rda")
load("Data/RFS_k_post.rda")

# 3.2 Define growth function that draws from above distributions -------------------------------------------------

# Stochastic growth function
FUN.growth <- function(z1, z, m.par) {
  linf <- sample(m.par[["Linf"]], 1)
  mu <- z + (linf - z) * (1 - exp(-sample(m.par[["k.mu"]], 1)))
  # mean size next year (z1)
  sig <- ifelse(mu > linf, 1, m.par[["k.sd"]])
  # sd about mean
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  # pdf for size z1
  return(p.den.grow)
}

# Constant growth function
FUN.growth.det <- function(z1, z, m.par) {
  linf <- mean(m.par[["Linf"]])
  mu <- z + (linf - z) * (1 - exp(-mean(m.par[["k.mu"]])))
  # mean size next year (z1)
  sig <- m.par[["k.sd"]]
  # sd about mean
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  # pdf for size z1
  return(p.den.grow)
}

# 3.3 Define parameters for the growth function ------------------------------------------------------------------
Growth.par <- list(
  # Mean Linf 90% credible interval
  Linf = sort(Linf.post), 
  # Standard deviation Linf 90% credible interval
  #Linf.sd = sort(Linf.sd.post)[-c(1:15,286:300)],
  # Mean growth rate 90% credible interval
  k.mu = sort(k.post),
  # Standard deviation of growth rate
  k.sd = 3)

# Make sure names match function call
names(Growth.par) <- c("Linf", "k.mu", "k.sd")


#-----------------------------------------------------------------------------------------------------------------
# 4.0 Adult survival kernel - Brooks 2020
#-----------------------------------------------------------------------------------------------------------------

# 4.1 Load survival posterior distributions ----------------------------------------------------------------------
load("Data/surv23_beta_pyear.RData")
load("Data/surv23_mu_pyear.RData")

# 4.2 Define survival function that draws from above distributions -------------------------------------------------

# Function for stochastic survival
FUN.surv <- function(z, m.par) {
  # linear predictor
  linear.p <- sample(m.par[["Surv.mu"]], 1) + sample(m.par[["Surv.beta"]], 1) * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

# Functions for constant survival at three levels
FUN.surv.90 <- function(z, m.par) {
  # linear predictor
  linear.p <- quantile(m.par[["Surv.mu"]], 0.90)[[1]] + quantile(m.par[["Surv.beta"]], 0.90)[[1]] * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

FUN.surv.70 <- function(z, m.par) {
  # linear predictor
  linear.p <- quantile(m.par[["Surv.mu"]], 0.70)[[1]] + quantile(m.par[["Surv.beta"]], 0.70)[[1]] * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

FUN.surv.50 <- function(z, m.par) {
  # linear predictor
  linear.p <- quantile(m.par[["Surv.mu"]], 0.50)[[1]] + quantile(m.par[["Surv.beta"]], 0.50)[[1]] * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

FUN.surv.30 <- function(z, m.par) {
  # linear predictor
  linear.p <- quantile(m.par[["Surv.mu"]], 0.30)[[1]] + quantile(m.par[["Surv.beta"]], 0.30)[[1]] * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

FUN.surv.10 <- function(z, m.par) {
  # linear predictor
  linear.p <- quantile(m.par[["Surv.mu"]], 0.10)[[1]] + quantile(m.par[["Surv.beta"]], 0.10)[[1]] * ((z - m.par[["SVL.mean"]])/m.par[["SVL.sd"]])
  # inv-logistic trans
  p <- 1 / (1 + exp(-linear.p))
  return(p)
}

# 4.3 Define parameters for the survival function ----------------------------------------------------------------
Surv.par <- list(
  ## Mean survival 90% credible interval
  sort(mu.post),
  ## Relationship with body size 90% credible interval
  sort(beta.post),
  # Center and scale SVL values from survival analysis
  SVL.mean = 62.6,
  SVL.sd = 9.73)

# Make sure names match function call
names(Surv.par) <- c("Surv.mu", "Surv.beta", "SVL.mean", "SVL.sd")


#-----------------------------------------------------------------------------------------------------------------
# 5.0 Combine kernels and parameters
#-----------------------------------------------------------------------------------------------------------------

# 5.1 Combine the survival and growth kernels --------------------------------------------------------------------

# Fully stochastic function
FUN.surv.z <- function(z1, z, m.par){
  return(FUN.surv(z, m.par) * FUN.growth(z1, z, m.par))
}

# Constant functions for three levels of survival (growth is always stochastic)
FUN.surv.z.90 <- function(z1, z, m.par){
  return(FUN.surv.90(z, m.par) * FUN.growth(z1, z, m.par))
}

FUN.surv.z.70 <- function(z1, z, m.par){
  return(FUN.surv.70(z, m.par) * FUN.growth(z1, z, m.par))
}

FUN.surv.z.50 <- function(z1, z, m.par){
  return(FUN.surv.50(z, m.par) * FUN.growth(z1, z, m.par))
}

FUN.surv.z.30 <- function(z1, z, m.par){
  return(FUN.surv.30(z, m.par) * FUN.growth(z1, z, m.par))
}

FUN.surv.z.10 <- function(z1, z, m.par){
  return(FUN.surv.10(z, m.par) * FUN.growth(z1, z, m.par))
}

# 5.2 Compile all parameter values needed to run the above functions ---------------------------------------------

# Combine all parameters from above
m.par <- c(Surv.par, Fec.par, Growth.par)


#-----------------------------------------------------------------------------------------------------------------
# 6.0 EXTRA FUNCTIONS
#-----------------------------------------------------------------------------------------------------------------

# 3.1 Combine above functions to create fecundity kernel ---------------------------------------------------------

# This is a function for fecundity that will pull either a short or long hydroperiod size function. However,
# I don't think it's easy to integrate into the model at this stage so opted to put this variation in the model
# at a later step.

FUN.fec.z <- function(z1, z, hydro, m.par){
  
  if(hydro == 0) {
    
    rep(0, length(z))
    
  } else if (hydro == 1) {
    
    (1/2) * FUN.prob.breed(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.short(z1, m.par)
    
  } else {
    
    (1/2) * FUN.prob.breed(z, m.par) * FUN.eggs(z, m.par) * FUN.prob.rec(1) * FUN.metas.long(z1, m.par)
    
  }
  
}

