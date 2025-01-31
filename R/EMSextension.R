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
# - 2. Use Bayesian approach to estimate posterior distribution of Z-matrix
# - 3. Iteratively ordinate a subset of the posterior Z-matrices
# - 4. For each ordinated matrix, store standardized scores for coherence and turnover 
#      to build a posterior distribution for each EMS statistic
# - 5. From all of the ordinated incidence matrices, generate a heat map of probable/true 
#      metacommunity structure. 

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
########## SIMULATION OF OCCUPANCY DATA ############
####################################################

# Establish basic parameters for the simulation:
N <- 10 # Number of species
K <- 50 # Number of sites
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


# Establish species-specific base-line occurrence
# I assume these base-line probabilities do not change over time.
b.spp <- rnorm(N, mu_b, sd_b)

# Establish species-specific detection probabilities
lrho.spp <- rnorm(N, mu_rho, sd_rho)
rho.spp <- AntiLogit(lrho.spp)

# Establish covariate effects per species:
# I will assume that occupancy is only potentially related to one covariate for each species.
betas.b <- rnorm(N, -1, 1)
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

#Check row and column sums:
zeros <- NULL
for(i in 1:nrow(Y)){
  if(sum(Y[i, ])==0) zeros <- c(zeros, i)
}
zeros
if(sum(zeros)>0) Y <- Y[-zeros, ]

zeros <- NULL
for(j in 1:ncol(Y)){
  if(sum(Y[, j])==0) zeros <- c(zeros, j)
}
zeros
if(sum(zeros)>0) Y <- Y[, -zeros]

# Visualize observed structure
library(metacom)
Y2 <- aperm(Y, c(2, 1)) # Now sites are rows and species are columns
Y2 <- ifelse(Y2 > 0, 1, 0)
ordY <- OrderMatrix(Y2)
colnames(ordY) <- as.factor(1:ncol(ordY))
quartz(height=6, width=3)
print(Matrix_Plot(ordY, xlab="Species", ylab="Sites"))


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
########## BAYESIAN OCCUPANCY MODELING #############
####################################################

jags_d <- list(X=X,
               Y=Y,
               K=ncol(Y),
               N=nrow(Y),
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

iter <- 30 # number of samples to draw from the posterior 

# Each chain has 300 observations (n.iter=3000, thinned by 10)
# Choose which iterations will be used from the posterior
samples <- sample(c(1:300), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(nrow(Y), ncol(Y), iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.out, Chain==chain & Iteration==samples[i])$value
  # Store this vector as a matrix
  z.mat <- NULL
  z.mat <- matrix(subset, nrow=nrow(Y), ncol=ncol(Y), byrow=F)
  # Add matrix to the array
  z.post[, , i] <- z.mat
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
######## CALCULATE METACOMMUNITY METRICS  ##########
####################################################

# Make sites rows and species columns:
z.post <- aperm(z.post, c(2, 1, 3))

# Ordinate all the z matrices:

library(metacom)
# Store the Ordinated Matrices:
z.ord <- array(0, dim=c(ncol(Y), nrow(Y), iter))
# Store the statistics for Coherence and Turnover.
Coher <- array(0, dim=c(iter, 5))
colnames(Coher) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
Turn <- array(0, dim=c(iter, 5))
colnames(Turn) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
Bound <- array(0, dim=c(iter, 3))
colnames(Bound) <- c("index", "pval", "df")

for(i in 1:iter){
  meta <- NULL
  meta <- metacommunity(z.post[, , i], method="r1", sims=1000, allow.empty=T)
  z.ord[, , i] <- meta[[1]]
  Coher[i, ] <- as.numeric(as.character(meta[[2]][1:5, ]))
  Turn[i, ] <- as.numeric(as.character(meta[[3]][1:5, ]))
  for(j in 1:3){
    Bound[i, j] <- meta[[4]][1,j]
  }
}

# Generate heat map of ordinated matrices:
quartz(height=6, width=3)
print(Matrix_HeatMap(z.ord, ylab="Sites"))

# Generate density plots of posterior metacommunity metrics:
Coher.plot <- ggplot(data.frame(Coher), aes(x=z))+
                geom_density(kernel="biweight")+
                labs(x="Coherence z-score", y="Density")+
                theme_classic()

Turn.plot <- ggplot(data.frame(Turn), aes(x=z))+
              geom_density(kernel="biweight")+
              labs(x="Turnover z-score", y="")+
              theme_classic()

Bound.plot <- ggplot(data.frame(Bound), aes(x=index))+
                geom_density(kernel="biweight")+
                labs(x="Morista's I", y="")+
                theme_classic()

library(gridExtra)
grob <- arrangeGrob(Coher.plot, Turn.plot, Bound.plot, nrow=1)
quartz(height=5, width=9)
print(grob)










#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
######################   UNUSED CODE BELOW   ###################### 
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------


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
  
  
  
  
  
  
  