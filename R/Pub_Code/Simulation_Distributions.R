
# Things to do:
# (1) Automate sequence of simulations:
# - alter distribution of covariate effects
# - nest alteration of covariate value distribution
# (2) Automate determination of metacommunity structure
# (3) Store results in a data.frame for easy access

# Permuations of covariate effects:
# fixed (0)
# Normal (0, 0.5)
# Normal (0, 1)
# Uniform (-1, 1)
# Uniform (-3, 3)

# Permutations of covariate values:
# Normal (0, 1)
# Normal (0, 5)
# Normal (2, 1)
# Normal (5, 2)
# Uniform (-1, 1)
# Uniform (-4, 4)
# Uniform (0, 3)
# Uniform (0, 10)

# Store normal and uniform parameters (means, sd; high, low)
cov.eff <- 
  c(0, 0.5,
    0, 1,
    -1, 1,
    -3, 3)
cov.val <-
  c(0, 1,
    0, 5, 
    2, 1,
    5, 2,
    -1, 1,
    -4, 4,
    0, 3,
    0, 10)

mat.eff <- matrix(cov.eff, ncol=2, byrow=T) # Columns are mean and sd or high and low
mat.val <- matrix(cov.val, ncol=2, byrow=T) 

# Establish basic parameters for the occupancy model:
N <- 12 # Number of species
K <- 75 # Number of sites
J <- 4  # Number of sampling replicates per site

# Establish some useful internal functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

# Base-line occurrence probability
b0 <- rep(Logit(0.6), N) # Fixed occurrence across species

# Establish species-specific detection probabilities
p0 <- rep(0.7, N) # Fixed detection prob.


#############################################
### SIMULATE COMMUNITY INCIDENCE MATRICES ###
#############################################

# Set number of iterations for each set of parameter values
iter <- 2
total.iter <- iter*nrow(mat.eff)*nrow(mat.val) # Total number of simulations

b.spp <- NULL
X <- NULL
  
z <- array(0, dim = c(N, K)) # This will hold occupancy states
Y <- array(0, dim = c(N, K, nrow(mat.eff), nrow(mat.val), iter)) # This will hold simulated observations
lpsi <- array(0, dim = c(N, K)) # This will hold the occupancy probabilities
psi <- array(0, dim = c(N, K))



for (it in 1:iter){
  for(i in 1:nrow(mat.eff)){
    # Draw species-specific covariate effects
    if(i < 3){ # Alter this number depending on the parameters in mat.eff
      b.spp <- rnorm(N, mat.eff[i, 1], mat.eff[i, 2])
    }else{
      b.spp <- runif(N, mat.eff[i, 1], mat.eff[i, 2])
    }
    
    for(j in 1:nrow(mat.val)){
      # Draw site-specific covariate values
      if(j < 5){ # Alter this number depending on the parameters in mat.val
        X <- rnorm(K, mat.val[j, 1], mat.val[j, 2])
      }else{
        X <- runif(K, mat.val[j, 1], mat.val[j, 2])
      }
      
      # Occupancy states
      for(n in 1:N){
        for(k in 1:K){
          lpsi[n, k] <- b0[n] + b.spp[n] * X[k]
          psi[n, k] <- AntiLogit(lpsi[n, k])
          
          z[n, k] <- rbinom(1, 1, psi[n, k])
          Y[n, k, i, j, it] <- rbinom(1, J, p0[n] * z[n, k])
        }
      }
    }
  }
}

####
# Now I need to figure out how to do the steps below (get rid of zero columns/rows), and then store in a list
# But how to get from array to list? 


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

# Calculate metacommunity metrics
library(metacom)
Y2 <- aperm(Y, c(2, 1)) # Now sites are rows and species are columns
Y2 <- ifelse(Y2 > 0, 1, 0)
Metacommunity(Y2)