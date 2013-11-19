# Author: JR Mihaljevic
# Date: November 2013

# Extending the "Elements of Metacommunity Structure" paradigm w/
# single time, multi-species occupancy models (multiple observations w/in one time unit)
# Goals: 
# - Generate a heat map of ordinated site-species incidence matrices
# - Generate a distribution of z-values for metacommunity statistics (e.g. coherence)

# This method will incorporate error in the probability of detection for each species into
# analysis of metacommunity structure. 

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# General Methods:
# - 1. Simulate data for a simple one-year, multi-species occupancy model
#       - assume that occupancy is structured along a single environmental covariate
# - 2. Use Bayesian approach to estimate posterior distribution of model parameters
# - 3. Build a simulation model to simulate many ordinated incidence matrices, drawing upon
#      the posterior distributions of the estimated parameters from the observed data-set
# - 4. For each ordinated matrix, store z-values of coherence and turnover to build a 
#      distribution of z-values for each EMS statistic
# - 5. From all of the ordinated incidence matrices, generate a heat map of probable 
#      metacommunity structure. 

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
########## SIMULATION OF OCCUPANCY DATA ############
####################################################

# Establish basic parameters for the simulation:
N <- 10 # Number of species
K <- 100 # Number of sites
J <- 4  # Number of sampling replicates per site

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}


# Establish community-level hyperparameters
p_b <- 0.5 # community-level probability of occurrence
mu_b <- Logit(p_b)
sd_b <- .75 # standard dev. b

p_rho <- 0.5 # community-level probability of detection
mu_rho <- Logit(p_rho)
sd_rho <- 1 


# Establish species-specific base-line random effects
# I assume these base-line effects do not change over time.
b.spp <- rnorm(N, mu_b, sd_b)

lrho.spp <- rnorm(N, mu_rho, sd_rho)
rho.spp <- AntiLogit(lrho.spp)

# Establish covariate effects per species:
# I will assume that occupancy is only potentially related to one covariate for each species.
betas.b <- rnorm(N, 0, 1)
#betas.b <- runif(N, -4, 1)

# Establish covariate values:
# I will assume these values do not change over the time-scale of the multiple observations.
X <- rnorm(K, 2, 3)


# Initial occupancy state
z <- array(0, dim = c(N, K)) # This will hold occupancy states
Y <- array(0, dim = c(N, K)) # This will hold simulated observations
lpsi <- array(0, dim = c(N, K)) # This will hold the occupancy probabilities
psi <- array(0, dim = c(N, K))

for(i in 1:N){
  for(k in 1:K){
    lpsi[i, k] <- b.spp[i] + betas.b[i] * X[k]
    psi[i, k] <- AntiLogit(lpsi[i, k])
    
    z[i, k] <- rbinom(1, 1, psi[i, k])
    Y[i, k] <- rbinom(1, J, rho.spp[i] * z[i, k])
  }
}

# Visualize the ordinated matrix of this simulated data:
Y2 <- aperm(Y, c(2, 1)) # Now sites are rows and species are columns
Y2 <- ifelse(Y2 > 0, 1, 0) # Make sure we just have zeros/ones
library(metacom)
#Check row and column sums:
zeros <- NULL
for(i in 1:nrow(Y2)){
  if(sum(Y2[i, ])==0) zeros <- c(zeros, i)
}
if(sum(zeros)>0) Y2 <- Y2[-zeros, ]

zeros <- NULL
for(j in 1:ncol(Y2)){
  if(sum(Y2[j, ])==0) zeros <- c(zeros, i)
}
if(sum(zeros)>0) Y2 <- Y2[, -zeros]

ordY <- OrderMatrix(Y2)

quartz(height=6, width=3)
print(Matrix_Plot(ordY, xlab="Species", ylab="Sites"))


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
########## BAYESIAN OCCUPANCY MODELING #############
####################################################

jags_d <- list(X=X,
               Y=Y,
               K=K,
               N=N,
               J=J)

# Set initial parameters:
# Z values (unobserved)
zinit <- Y
zinit <- ifelse(zinit > 0, 1, 0)

# Start the model
library(rjags)

params <- c("z") # This will give me species-level estimates for each of these

mod <- jags.model(file = "OccMod_SingleYear.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))

out <- coda.samples(mod, n.iter = 3000, variable.names = params, thin=10)
summary(out)
plot(out)

# Store output:
library(ggmcmc)
post.out <- ggs(out)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
####################################################

iter <- 40 # number of samples to draw from the posterior 

# Each chain has 300 observations (n.iter=3000, thinned by 10)
# Choose which iterations will be used from the posterior
samples <- sample(c(1:300), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(N, K, iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.out, Chain==chain & Iteration==samples[i])$value
  # Store this vector as a matrix
  z.mat <- NULL
  z.mat <- matrix(subset, nrow=N, ncol=K, byrow=F)
  # Add matrix to the array
  z.post[, , i] <- z.mat
}

# Make sites rows and species columns:
z.post <- aperm(z.post, c(2, 1, 3))

# Ordinate all the z matrices:
z.ord <- array(0, dim=c(K, N, iter))

library(metacom)
for(i in 1:iter){
  z.ord[, , i] <- OrderMatrix(z.post[, , i])
}

quartz(height=6, width=3)
print(Matrix_HeatMap(z.ord))


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
########## SIMULATE FROM BAYESIAN OUTPUT ###########
####################################################

# Now I'll use the posterior distrubtions of the estimated parameters to simulate new 
# incidence matrices

iter <- 50 # Number of iterations
Y.post <- array(0, dim = c(N, K, iter))
Y.ord <- array(0, dim = c(K, N, iter))
require(metacom)

for(i in 1:iter){
  
  b.spp.post <- NULL
  betas.b.post <- NULL
  rho.spp.post <- NULL
  
  for(n in 1:N){
    param <- paste("b0[",n,"]", sep="")
    b.spp.post[n] <- sample(subset(post.out, Parameter==param)$value, size=1)
    
    param <- paste("b[",n,"]", sep="")
    betas.b.post[n] <- sample(subset(post.out, Parameter==param)$value, size=1)
    
    param <- paste("p[",n,"]", sep="")
    rho.spp.post[n] <- sample(subset(post.out, Parameter==param)$value, size=1)
    
  }
  
  z.post <- array(0, dim = c(N, K)) # This will hold occupancy states
  lpsi.post <- array(0, dim = c(N, K)) # This will hold the occupancy probabilities
  psi.post <- array(0, dim = c(N, K))
  
  for(n in 1:N){
    for(k in 1:K){
      lpsi.post[n, k] <- b.spp.post[n] + betas.b.post[n] * X[k]
      psi.post[n, k] <- AntiLogit(lpsi.post[n, k])
      
      z.post[n, k] <- rbinom(1, 1, psi.post[n, k])
      Y.post[n, k, i] <- rbinom(1, J, rho.spp.post[n] * z.post[n, k])
    }
  }
}

Y.post <- ifelse(Y.post > 0, 1, 0)
Y.post <- aperm(Y.post, c(2, 1, 3))

# NEED TO FIGURE OUT A WAY TO DEAL WITH ROWS AND COLUMNS THAT SUM TO ZERO!!!
for(i in 1:iter){
  Y.ord[, , i] <- OrderMatrix(Y.post[, , i])
}

quartz(height=6, width=3)
print(Matrix_HeatMap(Y.ord, xlab="Species", ylab="Sites"))
  
  
  
  
  
  
  