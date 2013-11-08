# Simulating some data for multi-species, dynamic occupancy modeling.
# Author: JR Mihaljevic
# November 2013
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# I'm going to assume that one covariate affects each parameter:
# psi = prob occurrence (predicted by "b")
# gamma = prob colonization if not present in T-1 (predicted by "c")
# phi = prob persistence (predicted by "d")
# rho = prob detection (no covariates assumed here)

# Establish basic parameters for the simulation:
n <- 6 # Number of species
K <- 100 # Number of sites
J <- 5  # Number of sampling replicates per site, per time period
Ts <- 4 # Number of sampling time points
ncov <- 3 # Number of covariates in model

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Establish community-level hyperparameters
p_b <- 0.7 # community-level probability of occurrence
mu_b <- Logit(p_b)
sd_b <- 1 # standard dev. b

p_c <- 0.3 # community-level probability of colonization
mu_c <- Logit(p_c)
sd_c <- 2

p_d <- 0.6 # community-level probability of persistence per year
mu_d <- Logit(p_d)
sd_d <- 2

p_rho <- 0.8 # community-level probability of detection
mu_rho <- Logit(p_rho)
sd_rho <- 1 

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Establish species-specific base-line random effects
# I assume these base-line effects do not change over time.
b.spp <- rnorm(n, mu_b, sd_b)
c.spp <- rnorm(n, mu_c, sd_c)
d.spp <- rnorm(n, mu_d, sd_d)
lrho.spp <- rnorm(n, mu_rho, sd_rho)
rho.spp <- AntiLogit(lrho.spp)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Establish covariate effects per species:
# I will assume that each parameter is only related to one of the three covariates.
betas.b <- mat.or.vec(nr=n, nc=ncov)
for(i in 1:n){
  for(j in 1:ncov){
    if(j==1) {betas.b[i, j] <- rnorm(1, 0, 2)}
    else{betas.b[i, j] <- rnorm(1, 0, 0.1)}
  }
}

betas.c <- mat.or.vec(nr=n, nc=ncov)
for(i in 1:n){
  for(j in 1:ncov){
    if(j==2) {betas.c[i, j] <- rnorm(1, 0, 2)}
    else{betas.c[i, j] <- rnorm(1, 0, 0.1)}
  }
}

betas.d <- mat.or.vec(nr=n, nc=ncov)
for(i in 1:n){
  for(j in 1:ncov){
    if(j==3) {betas.d[i, j] <- rnorm(1, 0, 2)}
    else{betas.d[i, j] <- rnorm(1, 0, 0.1)}
  }
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Establish covariate values, changing through time:
X <- array(0, dim=c(K, ncov, Ts))

for(k in 1:K){
  for(cov in 1:ncov){
    for(t in 1:Ts){
      X[k, cov, t] <- rnorm(1, 0, 10)
    }
  }
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Initial occupancy state
z <- array(0, dim = c(n, K, Ts)) # This will hold occupancy states
Y <- array(0, dim = c(n, K, Ts)) # This will hold simulated observations
lpsi <- array(0, dim = c(n, K, Ts)) # This will hold the occupancy probabilities
psi <- array(0, dim = c(n, K, Ts))

for(i in 1:n){
  for(k in 1:K){
    lpsi[i, k, 1] <- b.spp[i] + betas.b[i, ] %*% X[k, ,1]
    psi[i, k, 1] <- AntiLogit(lpsi[i, k, 1])
    z[i, k, 1] <- rbinom(1, 1, psi[i, k, 1])
    Y[i, k, 1] <- rbinom(1, J, rho.spp[i] * z[i, k, 1])
  }
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Subsequent occupancy states:
lphi <- array(0, dim = c(n, K, Ts-1)) # This will hold the persistence probabilities
phi <- array(0, dim = c(n, K, Ts))
lgam <- array(0, dim = c(n, K, Ts-1)) # This will hold the colonization probabilities
gam <- array(0, dim = c(n, K, Ts-1))

for(i in 1:n){
  for(t in 1:(Ts-1)){
    for (k in 1:K) { # K=number of sites
      lgam[i,k,t] <- c.spp[i] + betas.c[i, ] %*% X[k, ,t]
      gam[i,k,t] <-  AntiLogit(lgam[i,k,t])
      
      lphi[i,k,t] <- d.spp[i] + betas.d[i, ] %*% X[k, ,t]
      phi[i,k,t] <-  AntiLogit(lphi[i,k,t])
      
      # Psi[t+1] = persist * occupy[t] + colonize*(1-occupy[t])
      psi[i,k,t+1] <- phi[i,k,t] * psi[i,k,t] + gam[i,k,t] * (1 - psi[i,k,t])
      
      z[i,k,t+1] <- rbinom(1, 1, phi[i,k,t] * z[i,k,t] + gam[i,k,t] * (1 - z[i,k,t]))
      
      Y[i,k,t+1] <- rbinom(1, J, rho.spp[i] * z[i, k, t+1])
    } 
  }
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Enter data into list:

jags_d <- list(x=X,
               y=Y,
               K=K,
               ncovs=ncov,
               T=Ts,
               n=n,
               J=J)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Set initial parameters:

# Z values (unobserved):
zinit <- Y
zinit <- ifelse(zinit > 0, 1, 0)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Start the model
library(rjags)

mod <- jags.model(file = "OccupancyModel_Draft.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))

out <- coda.samples(mod, n.iter = 5000, variable.names = c("b0", "b"))
summary(out)
plot(out)





