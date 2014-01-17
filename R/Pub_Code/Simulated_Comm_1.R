# Simulate occupancy data for:
# Figure 1
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


# Establish community-level hyperparameters
p_b <- 0.5 # community-level probability of occurrence
mu_b <- Logit(p_b)
sd_b <- 0.5 # standard dev. b

p_rho <- 0.6 # community-level probability of detection
mu_rho <- Logit(p_rho)
sd_rho <- 0.75 


# Establish species-specific base-line occurrence
# I assume these base-line probabilities do not change over time.
#b.spp <- rnorm(N, mu_b, sd_b)
b.spp <- rep(0.4, N)

# Establish species-specific detection probabilities
lrho.spp <- rnorm(N, mu_rho, sd_rho)
rho.spp <- AntiLogit(lrho.spp)

# Establish covariate effects per species:
# I will assume that occupancy is only potentially related to one covariate for each species.
# betas.b <- rnorm(N, 0, 0.75)
betas.b <- runif(N, -1.5, 1.5)

# Establish covariate values:
# I will assume these values do not change over the time-scale of the multiple observations.
# X <- rnorm(K, 0, 1.5)
X <- runif(K, -4, 4)


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

