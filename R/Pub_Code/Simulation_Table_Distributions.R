# Create a table to show how changes in...
# - Baseline occ. prob
# - Distribution of covariate effects
# - Distribution of covariate values
# ... affect the metacommunity structure of a simulated metacommunity. 
# Author: JR Mihaljevic
# January 2014


####################################################
########## SIMULATION OF OCCUPANCY DATA ############
####################################################

# Establish basic parameters for the simulation:
N <- 12 # Number of species
K <- 75 # Number of sites
J <- 4  # Number of sampling replicates per site

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}


# Establish community-level hyperparameters (if using normal distribution)
# community-level probability of occurrence
# p_b <- 0.3 
# mu_b <- Logit(p_b)
# sd_b <- 0.75 # standard dev. b


# Establish species-specific base-line occurrence
# I assume these base-line probabilities do not change over time.
b.spp <- rep(Logit(0.6), N) # Fixed occurrence
# b.spp <- rnorm(N, mu_b, sd_b) # Normal distribution 
# b.spp <- runif(N, 0, 1) # Uniform distribution

# # Establish species-specific detection probabilities
rho.spp <- rep(0.7, N) # Fixed detection prob. 


# Establish covariate effects per species:
# I will assume that occupancy is only potentially related to one covariate for each species.
# betas.b <- rep(0, N)
# betas.b <- rnorm(N, 0, .5)
 betas.b <- runif(N, -1, 1)

# For nested patterns:
#betas.b <- runif(N, -1, 1)

# Establish covariate values:
# I will assume these values do not change over the time-scale of the multiple observations.
#X <- rnorm(K, 0, 2) # Normal dist.
X <- runif(K, -3, 3) # Uniform dist.

# For nested patterns:
#X <- runif(K, 0, 3)
#X <- rnorm(K, 2, 1)

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

# Calculate metacommunity metrics
Metacommunity(Y2)
