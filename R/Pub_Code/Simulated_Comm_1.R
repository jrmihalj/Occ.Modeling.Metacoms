# Determining bias in assigning metacommunity structure due to detection error
# Does an occupancy model do a better job at assigning 'known' structure of simulated communities?
# Author: JR Mihaljevic
# January 2014

##########################################
##########       TO DO:       ############
##########################################

# 1. CHANGE THE NUMBER OF ITERATIONS (HOW MANY UNIQUE COMMUNITIES TO SIMULATE)
#     - DEFAULT AT 1000

# 2. CHANGE THE DETECTION PROBABILITIES
#     - DEFAULT AT 70%, FIXED AMONG ALL SPECIES

# 3. CHANGE THE NUMBER OF Z_POSTERIOR MATRICES TO USE
#     - DEFAULT AT 1000 (NOT EASY TO CHANGE AT THIS POINT)

# 4. SET YOUR WORKING DIRECTORY TO WHERE THE JAGS OCCUPANCY MODEL IS LOCATED

##################################################
########## ESTABLISH USEFUL FUNCTIONS ############
##################################################

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

# Function to calculate Rhat for MCMC chain convergence
# Adapted from the ggs_Rhat() function in the 'ggmcmc' package
get_Rhat <- function (D, family = NA, scaling = 1.5) 
{
  if (attributes(D)$nChains < 2) {
    stop("At least two chains are required")
  }
  if (!is.na(family)) {
    D <- get_family(D, family = family)
  }
  psi.dot <- ddply(D, .(Parameter, Chain), summarize, psi.dot = mean(value), 
                   .parallel = attributes(D)$parallel)
  psi.j <- ddply(D, .(Parameter), summarize, psi.j = mean(value), 
                 .parallel = attributes(D)$parallel)
  b.df <- merge(psi.dot, psi.j)
  attr(b.df, "nIterations") <- attributes(D)$nIterations
  b.df <- cbind(b.df, nIterations = attributes(D)$nIterations)
  B <- ddply(b.df, .(Parameter), summarize, B = var(psi.j - 
                                                      psi.dot) * nIterations, .parallel = attributes(D)$parallel)
  s2j <- ddply(D, .(Parameter, Chain), summarize, s2j = var(value), 
               .parallel = attributes(D)$parallel)
  W <- ddply(s2j, .(Parameter), summarize, W = mean(s2j), .parallel = attributes(D)$parallel)
  BW <- merge(B, W)
  BW <- cbind(BW, nIterations = attributes(D)$nIterations)
  BW <- ddply(BW, .(Parameter), transform, wa = ((((nIterations - 
                                                      1)/nIterations) * W) + ((1/nIterations) * B)), .parallel = attributes(D)$parallel)
  BW <- ddply(BW, .(Parameter), transform, Rhat = sqrt(wa/W), 
              .parallel = attributes(D)$parallel)
  return(unique(BW$Rhat))
}

##############################################
########## LOAD REQUIRED PACKAGES ############
##############################################

require(metacom)
require(rjags)
require(ggmcmc)

####################################################
########## SIMULATION OF OCCUPANCY DATA ############
####################################################
# Set number of iterations (How many unique Z and Y matrices to create?)
iter <- 1000

# Establish basic parameters for the simulation:
N <- 12 # Number of species
K <- 75 # Number of sites
J <- 3  # Number of sampling replicates per site

# Base-line occurrence probability
# Fixed occurrence across species, in logit space
b0 <- rep(Logit(0.6), N) # Change if desired

# Community-level detection probability (hyper-parameters)
# p_mean <- 0.7
# mu_p <- Logit(p_mean) # Logit space
# sd_p <- 0.75 # standard dev.

# Fixed detection across species, in logit space
p0 <- rep(0.7, N) # Change if desired

# Storage
z <- array(0, dim = c(N, K, iter)) # 'True' occupancy states
Y <- array(0, dim = c(N, K, iter)) # Simulated observations
X <- array(0, dim = c(K, iter))    # Covariate values

for(i in 1:iter){

  b.spp <- rnorm(N, 0, 0.5)
  X[, i] <- rnorm(K, 0, 5)
  
#   lp0 <- rnorm(N, mu_p, sd_p) # Use if detection assigned by hyperparameters
#   p0 <- AntiLogit(lp0)
      
  #### Occupancy states ####
  # Storage:
  lpsi <- array(0, dim = c(N, K)) # occupancy probabilities (see below)
  psi <- array(0, dim = c(N, K))
  for(n in 1:N){
    for(k in 1:K){
      lpsi[n, k] <- b0[n] + b.spp[n] * X[k, i] # Covariate effects on occurrence 
      psi[n, k] <- AntiLogit(lpsi[n, k])
          
      z[n, k, i] <- rbinom(1, 1, psi[n, k]) # True Occupancy
      Y[n, k, i] <- rbinom(1, J, p0[n] * z[n, k, i]) # Observed Occupancy
      }
    }
}

#############################################
####  STORE SEPARATE MATRICES IN A LIST  ####
#############################################

# Storage:
Z_mats <- list()
Yobs_mats <- list()

for (i in 1:iter){
  Z_mats[[i]] <- z[, , i]
  Yobs_mats[[i]] <- Y[, , i]
}


#############################################
### CHECK EACH MATRIX FOR ROW/COLUMN SUMS ###
#############################################

# Check for rows and columns that sum to zero.
# These are not allowed in the 'metacom' package or 
# the associated functions in the 'vegan' package.

for(i in 1:iter){
  # Check column sums
  if(any(colSums(Z_mats[[i]]) == 0) == TRUE){
    Z_mats[[i]] <- Z_mats[[i]][, -which(colSums(Z_mats[[i]]) == 0)]
  }
  if(any(colSums(Yobs_mats[[i]]) == 0) == TRUE){
    Yobs_mats[[i]] <- Yobs_mats[[i]][, -which(colSums(Yobs_mats[[i]]) == 0)]
  }
  # Check row sums
  if(any(rowSums(Z_mats[[i]]) == 0) == TRUE){
    Z_mats[[i]] <- All_mats[[i]][-which(rowSums(Z_mats[[i]]) == 0), ]
  }
  if(any(rowSums(Yobs_mats[[i]]) == 0) == TRUE){
    Yobs_mats[[i]] <- Yobs_mats[[i]][-which(rowSums(Yobs_mats[[i]]) == 0), ]
  }
}

#############################################
# CALCULATE AND STORE METACOMMUNITY METRICS #
##########    FOR EACH MATRIX    ############
#############################################

#Storage of metacommunity metrics
Coher <- array(0, dim=c(iter, 5, 2)) # 2: One for Z_mats, one for Yobs_mats
colnames(Coher) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
Turn <- array(0, dim=c(iter, 5, 2))
colnames(Turn) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
Bound <- array(0, dim=c(iter, 3, 2))
colnames(Bound) <- c("index", "pval", "df")

for(i in 1:iter){
  # Need to transpose and convert to pres/abs
  mat_Z <- NULL
  mat_Y <- NULL
  
  mat_Z <- aperm(Z_mats[[i]], c(2,1))
  mat_Y <- aperm(Yobs_mats[[i]], c(2,1))
  
  mat_z <- ifelse(mat_Z > 0, 1, 0)
  mat_Y <- ifelse(mat_Y > 0, 1, 0)
  
  # Store all metacommunity metrics
  meta_Z <- NULL
  meta_Y <- NULL
  
  meta_Z <- Metacommunity(mat_Z, method="r1", sims=1000)
  meta_Y <- Metacommunity(mat_Y, method="r1", sims=1000)
  
  # Store Coherence Metrics
  Coher[i, ,1] <- as.numeric(as.character(meta_Z[[2]][1:5, ]))
  Coher[i, ,2] <- as.numeric(as.character(meta_Y[[2]][1:5, ]))
  # Store Turnover Metrics
  Turn[i, ,1] <- as.numeric(as.character(meta_Z[[3]][1:5, ]))
  Turn[i, ,2] <- as.numeric(as.character(meta_Y[[3]][1:5, ]))
  # Store Boundary Clumping Metrics
  for(j in 1:3){
    Bound[i, j, 1] <- meta_Z[[4]][1,j]
    Bound[i, j, 2] <- meta_Y[[4]][1,j]
  }
}

##############################################
###   DETERMINE METACOMMUNITY STRUCTURE   ####
###########    FOR EACH MATRIX    ############
##############################################

# Storage:
Structure <- array(0, dim=c(iter, 2))

for(i in 1:iter){
  for(j in 1:2){
    if(Coher[i, 3, j] > 0.05){
      Structure[i, j] <- "Random"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] > Coher[i, 4, j]){
      Structure[i, j] <- "Checkerboard"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] > 0.05 & Turn[i, 1, j] < Turn[i, 4, j]){
      Structure[i, j] <- "Quasi-Nested"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] > 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 1, j] > 1 & Bound[i, 2, j] <= 0.05){
      Structure[i, j] <- "Quasi-Clementsian"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] > 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 1, j] < 1 & Bound[i, 2, j] <= 0.05){
      Structure[i, j] <- "Quasi-EvenSpaced"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] > 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 2, j] > 0.05){
      Structure[i, j] <- "Quasi-Gleasonian"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] <= 0.05 & Turn[i, 1, j] < Turn[i, 4, j]){
      Structure[i, j] <- "Nested"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] <= 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 1, j] > 1 & Bound[i, 2, j] <= 0.05){
      Structure[i, j] <- "Clementsian"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] <= 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 1, j] < 1 & Bound[i, 2, j] <= 0.05){
      Structure[i, j] <- "EvenSpaced"
    }
    
    if(Coher[i, 3, j] <= 0.05 & Coher[i, 1, j] < Coher[i, 4, j] & 
         Turn[i, 3, j] <= 0.05 & Turn[i, 1, j] > Turn[i, 4, j] &
         Bound[i, 2, j] > 0.05){
      Structure[i, j] <- "Gleasonian"
    }
  }
}

which(Structure[, 1] != Structure[, 2])

# Determine the proportion of instances where Yobs does not match Z
length(which(Structure[, 1] != Structure[, 2])) / iter

############################################
######## BAYESIAN ESTIMATES FOR Z ##########
############################################

# Storage:
post.z <- list()

for(i in 1:iter){
  
  jags_d <- list(X=X[, i],
                 Y=Yobs_mats[[i]],
                 K=ncol(Yobs_mats[[i]]),
                 N=nrow(Yobs_mats[[i]]),
                 J=J)
  
  # Set initial parameters:
  # Z values (unobserved)
  
  zinit <- ifelse(Yobs_mats[[i]] > 0, 1, 0)
  
  # Start the model
  
  params <- c("lpsiMean", "lpsiSD", "z")
  
  mod <- jags.model(file = "OccMod_SingleYear.txt", 
                    data = jags_d, n.chains = 3, n.adapt=1000,
                    inits = list(z=zinit))
  update(mod, n.iter=5000) # 5000 burn-in
  
  out <- coda.samples(mod, n.iter = 10000, variable.names = params, thin=10)
  
  # Check for convergence for these two parameters:
  post.psimean <- ggs(out, family="lpsiMean")
  post.psisd <- ggs(out, family="lpsiSD")
  Rhat_mean <- get_Rhat(post.psimean)
  Rhat_sd <- get_Rhat(post.psisd)
  
  # If not converged, update model until convergence is achieved
  if(Rhat_mean > 1.1 | Rhat_sd > 1.1){
    repeat{
      update(mod, n.iter=5000) # An extra 5,000 burn-in
      
      out <- coda.samples(mod, n.iter = 10000, variable.names = params, thin=10)
      
      post.psimean <- ggs(out, family="lpsiMean")
      post.psisd <- ggs(out, family="lpsiSD")
      Rhat_mean <- get_Rhat(post.psimean)
      Rhat_sd <- get_Rhat(post.psisd)
      
      if(Rhat_mean < 1.1 & Rhat_sd < 1.1) {break}
    }
  }
  
  #Store z_ij values. 
  post.z[[i]] <- ggs(out, family="z")
  
}


####################################################
######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
####################################################

# Create a list with 'iter' elements.
# Each element is another list with the 1000 Z posterior matrices for each iter
# I randomly choose which chain to use for each posterior estimate

all_Zposts <- list()
for(i in 1:iter){
  #Storage:
  Zpost_mats <- list()
  
  for(j in 1:1000){
    # Randomly choose from which chain the sample will originate
    chain <- NULL
    chain <- sample(c(1:3), 1)
    
    # Create a subset vector of the values for all z[n, k]:
    subset <- NULL
    subset <- subset(post.z[[i]], Chain==chain & Iteration==paste(j))$value
    # Store this vector as a matrix
    z.mat <- NULL
    z.mat <- matrix(subset, nrow=nrow(Yobs_mats[[i]]), ncol=ncol(Yobs_mats[[i]]), byrow=F)
    
    # Add matrix to the array
    Zpost_mats[[j]] <- z.mat
    
    # Check for zero sum columns/rows
    if(any(colSums(Zpost_mats[[j]]) == 0) == TRUE){
      Zpost_mats[[j]] <- Zpost_mats[[j]][, -which(colSums(Zpost_mats[[j]]) == 0)]
    }
    if(any(rowSums(Zpost_mats[[j]]) == 0) == TRUE){
      Zpost_mats[[j]] <- Zpost_mats[[j]][-which(rowSums(Zpost_mats[[j]]) == 0), ]
    }
    
  }
  # Store the Zpost_mats list within the overall list
  all_Zposts[[i]] <- Zpost_mats
}

#############################################
# CALCULATE AND STORE METACOMMUNITY METRICS #
##########    FOR EACH Z_POST    ############
#############################################

#Storage of metacommunity metrics
Coher.Zpost <- array(0, dim=c(iter, 5, 1000))
colnames(Coher.Zpost) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
Turn.Zpost <- array(0, dim=c(iter, 5, 1000))
colnames(Turn.Zpost) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
Bound.Zpost <- array(0, dim=c(iter, 3, 1000))
colnames(Bound.Zpost) <- c("index", "pval", "df")

for(i in 1:iter){
  for(j in 1:1000){
    mat <- NULL
    mat <- all_Zposts[[i]][[j]] # Extract relevant matrix
    mat <- aperm(mat, c(2,1)) # Transpose
    
    # Store all metacommunity metrics
    meta <- NULL
    
    meta <- Metacommunity(mat, method="r1", sims=1000)
    
    # Store Coherence Metrics
    Coher.Zpost[i, ,j] <- as.numeric(as.character(meta_Z[[2]][1:5, ]))
    
    # Store Turnover Metrics
    Turn.Zpost[i, ,j] <- as.numeric(as.character(meta_Z[[3]][1:5, ]))
    
    # Store Boundary Clumping Metrics
    for(k in 1:3){
      Bound.Zpost[i, k, j] <- meta_Z[[4]][1,k]
    }
    
  }  
}

##############################################
###   DETERMINE METACOMMUNITY STRUCTURE   ####
###########    FOR EACH Z_POST    ############
##############################################

# Store output:
Structure.Zpost <- array(0, dim=c(iter, 1000))

for(i in 1:iter){
  for(j in 1:1000){
    if(Coher.Zpost[i, 3, j] > 0.05){
      Structure.Zpost[i, j] <- "Random"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] > Coher.Zpost[i, 4, j]){
      Structure.Zpost[i, j] <- "Checkerboard"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] > 0.05 & Turn.Zpost[i, 1, j] < Turn.Zpost[i, 4, j]){
      Structure.Zpost[i, j] <- "Quasi-Nested"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] > 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 1, j] > 1 & Bound.Zpost[i, 2, j] <= 0.05){
      Structure.Zpost[i, j] <- "Quasi-Clementsian"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] > 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 1, j] < 1 & Bound.Zpost[i, 2, j] <= 0.05){
      Structure.Zpost[i, j] <- "Quasi-EvenSpaced"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] > 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 2, j] > 0.05){
      Structure.Zpost[i, j] <- "Quasi-Gleasonian"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] <= 0.05 & Turn.Zpost[i, 1, j] < Turn.Zpost[i, 4, j]){
      Structure.Zpost[i, j] <- "Nested"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] <= 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 1, j] > 1 & Bound.Zpost[i, 2, j] <= 0.05){
      Structure.Zpost[i, j] <- "Clementsian"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] <= 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 1, j] < 1 & Bound.Zpost[i, 2, j] <= 0.05){
      Structure.Zpost[i, j] <- "EvenSpaced"
    }
    
    if(Coher.Zpost[i, 3, j] <= 0.05 & Coher.Zpost[i, 1, j] < Coher.Zpost[i, 4, j] & 
         Turn.Zpost[i, 3, j] <= 0.05 & Turn.Zpost[i, 1, j] > Turn.Zpost[i, 4, j] &
         Bound.Zpost[i, 2, j] > 0.05){
      Structure.Zpost[i, j] <- "Gleasonian"
    }
  }
}

##########################################
###   DETERMINE NUMBER OF DEVIATIONS   ###
########    FOR EACH POSTERIOR    ########
##########################################

# Storage:
percent.dev <- vector()
for(i in 1:iter){
  percent.dev[i] <- length(which(Structure.Zpost[i, ] != Structure[i, 1])) / 1000
}

summary(percent.dev)