# Simulating different metacommunity structures using an occupancy model framework

# Alters the distribution of species-specific covariate effects and
# the distributions of covariate values across sites:

# Permuations of covariate effects:
# Normal (0, 0.5)
# Normal (0, 1)
# Uniform (-3, 0)

# Permutations of covariate values (for each covariate effect):
# Normal (0, 1)
# Normal (0, 5)
# Normal (2, 1)
# Normal (5, 2)
# Uniform (-1, 1)
# Uniform (-4, 4)
# Uniform (0, 3)
# Uniform (0, 10)

# Author: JR Mihaljevic
# Date: March 2014


###########################
####   LOAD PACKAGES   ####
###########################


require(metacom)
require(reshape2)

################################
####   REQUIRED FUNCTIONS   ####
################################

# Establish some useful internal functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

###############################
####   SIMULATION SET-UP   ####
###############################

# Store normal and uniform parameters (means, sd; high, low)
cov.eff <- 
  c(0, 0.5,
    0, 1,
    -3, 0)
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

# Base-line occurrence probability
b0 <- rep(Logit(0.6), N) # Fixed occurrence across species, in logit space

# Establish species-specific detection probabilities
p0 <- rep(0.7, N) # Fixed detection prob.


############################
####   RUN SIMULATION   ####
############################

# Set number of iterations for each set of parameter values
iter <- 1

# Storage for Structure identities
Structure <- array("A", dim=c(nrow(mat.eff), nrow(mat.val), iter))


for (i in 1:nrow(mat.eff)){
  for(j in 1:nrow(mat.val)){
    for(it in 1:iter){
      
      #############################################
      ### SIMULATE COMMUNITY INCIDENCE MATRICES ###
      #############################################
      
      b.spp <- NULL
      X <- NULL
      
      # Draw species-specific covariate effects
      if(i < 3){ # Alter this number depending on the parameters in mat.eff
        b.spp <- rnorm(N, mat.eff[i, 1], mat.eff[i, 2])
      }else{
        b.spp <- runif(N, mat.eff[i, 1], mat.eff[i, 2])
      }
      
      # Draw site-specific covariate values
      if(j < 5){ # Alter this number depending on the parameters in mat.val
        X <- rnorm(K, mat.val[j, 1], mat.val[j, 2])
      }else{
        X <- runif(K, mat.val[j, 1], mat.val[j, 2])
      }
      
      # Set storage arrays:
      Z <- array(0, dim = c(N, K)) # 'True' occupancy states
      Y <- array(0, dim = c(N, K)) # Simulated observations
      lpsi <- array(0, dim = c(N, K)) # occupancy probabilities (see below)
      psi <- array(0, dim = c(N, K))
      
      # Occupancy states
      for(n in 1:N){
        for(k in 1:K){
          lpsi[n, k] <- b0[n] + b.spp[n] * X[k]
          psi[n, k] <- AntiLogit(lpsi[n, k])
          
          Z[n, k] <- rbinom(1, 1, psi[n, k])
          Y[n, k] <- rbinom(1, J, p0[n] * Z[n, k])
        }
      }
    
      ########################################
      ### CHECK MATRIX FOR ROW/COLUMN SUMS ###
      ########################################
      
      # Check column sums
      if(any(colSums(Z) == 0) == TRUE){
        Z <- Z[, -which(colSums(Z) == 0)]
      }
      
      # Check row sums
      if(any(rowSums(Z) == 0) == TRUE){
        Z <- Z[-which(rowSums(Z) == 0), ]
      }
      
      #############################################
      # CALCULATE AND STORE METACOMMUNITY METRICS #
      #############################################
      
      meta <- NULL
      Coher <- NULL
      Turn <- NULL
      Bound <- NULL
      
      meta <- tryCatch(Metacommunity(Z, method="r1", sims=1000),
                       error = function(e){"ERROR"})
      
      if(meta == "ERROR"){
        Structure[i, j, it] <- "ERROR"
      }else {
        # Store metacommunity metrics
        Coher <- as.numeric(as.character(meta[[2]][1:5, ]))
        Turn <- as.numeric(as.character(meta[[3]][1:5, ]))
        for(p in 1:3){
          Bound[p] <- meta[[4]][1,p]
        }
        
        ##############################################
        ###   DETERMINE METACOMMUNITY STRUCTURE   ####
        ##############################################
        
        if(Coher[3] > 0.05){
          Structure[i, j, it] <- "Random"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] > Coher[4]){
          Structure[i, j, it] <- "Checkerboard"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] > 0.05 & Turn[1] < Turn[4]){
          Structure[i, j, it] <- "Quasi-Nested"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] > 0.05 & Turn[1] > Turn[4] &
             Bound[1] > 1 & Bound[2] <= 0.05){
          Structure[i, j, it] <- "Quasi-Clementsian"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] > 0.05 & Turn[1] > Turn[4] &
             Bound[1] < 1 & Bound[2] <= 0.05){
          Structure[i, j, it] <- "Quasi-EvenSpaced"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] > 0.05 & Turn[1] > Turn[4] &
             Bound[2] > 0.05){
          Structure[i, j, it] <- "Quasi-Gleasonian"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] <= 0.05 & Turn[1] < Turn[4]){
          Structure[i, j, it] <- "Nested"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] <= 0.05 & Turn[1] > Turn[4] &
             Bound[1] > 1 & Bound[2] <= 0.05){
          Structure[i, j, it] <- "Clementsian"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] <= 0.05 & Turn[1] > Turn[4] &
             Bound[1] < 1 & Bound[2] <= 0.05){
          Structure[i, j, it] <- "EvenSpaced"
        }
        
        if(Coher[3] <= 0.05 & Coher[1] < Coher[4] & 
             Turn[3] <= 0.05 & Turn[1] > Turn[4] &
             Bound[2] > 0.05){
          Structure[i, j, it] <- "Gleasonian"
        }
        
      }
    
    }
  }
}

Structure.melted <- melt(Structure) # Now these can be rowbinded





