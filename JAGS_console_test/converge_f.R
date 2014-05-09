# Testing convergence of JAGS model for a variety 
# of metacommunity structures/iterations

##########################################
##########       TO DO:       ############
##########################################

# 1. CHANGE THE NUMBER OF ITERATIONS

# 2. SET YOUR WORKING DIRECTORY TO WHERE THE JAGS OCCUPANCY MODEL IS LOCATED
#    - This is where the JAGS output will be dumped as well

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

require(metacom) # and all dependencies
require(coda) # and all dependencies
require(ggmcmc) # and all dependencies

#######################
#  SET UP SIMULATION  #
#######################

# Set number of iterations (How many unique Z and Y matrices to create?)
iter <- 2

# Fixed detection probability across species.
# The simulation will cycle through each of these:
p0_means <- c(0.50)
sd_p <- 0.75 #St. Dev. fixed for all simulations

# Storage (Outputs):
ERROR <- NULL
Converge <- NULL
##########################
###  START SIMULATION  ###
##########################

for(i in 1:iter){
  # Establish basic parameters for the simulation:
  N <- 12 # Number of species
  K <- 75 # Number of sites
  J <- 4  # Number of sampling replicates per site
  
  mu_p <- NULL
  mu_p <- Logit(p0_means)
  
  # Base-line occurrence probability
  b0 <- rep(Logit(0.6), N) # Change if desired
  
  b.spp <- NULL
  X <- NULL
  lp0 <- NULL
  p0 <- NULL
  
  b.spp <- rnorm(N, 0, runif(1, 0.5, 1)) # mean held at 0
  X <- rnorm(K, runif(1, 0, 2), runif(1, 2, 5))
  lp0 <- rnorm(N, mu_p, sd_p) # Use if detection assigned by hyperparameters
  p0 <- AntiLogit(lp0)
  
  #### Occupancy states ####
  # Storage:
  Z <- array(0, dim = c(N, K)) # 'True' occupancy states
  Y <- array(0, dim = c(N, K)) # Simulated observations
  lpsi <- array(0, dim = c(N, K)) # occupancy probabilities (see below)
  psi <- array(0, dim = c(N, K))
  
  #### Simulate occupancy states ####
  for(n in 1:N){
    for(k in 1:K){
      lpsi[n, k] <- b0[n] + b.spp[n] * X[k] # Covariate effects on occurrence 
      psi[n, k] <- AntiLogit(lpsi[n, k])
      
      Z[n, k] <- rbinom(1, 1, psi[n, k]) # True Occupancy
      Y[n, k] <- rbinom(1, J, p0[n] * Z[n, k]) # Observed Occupancy
    }
  }
  
  ########################################
  ### CHECK MATRIX FOR ROW/COLUMN SUMS ###
  ########################################
  
  # Check column sums
  if(any(colSums(Z) == 0) == TRUE){
    Z <- Z[, -which(colSums(Z) == 0)]
  }
  if(any(colSums(Y) == 0) == TRUE){
    Y <- Y[, -which(colSums(Y) == 0)]
  }
  # Check row sums
  if(any(rowSums(Z) == 0) == TRUE){
    Z <- Z[-which(rowSums(Z) == 0), ]
  }
  if(any(rowSums(Y) == 0) == TRUE){
    Y <- Y[-which(rowSums(Y) == 0), ]
  }
  #############################################
  # CALCULATE AND STORE METACOMMUNITY METRICS #
  ##########    FOR EACH MATRIX    ############
  #############################################
  
  mat_Z <- NULL
  mat_Y <- NULL
  # Transpose to fit requirements of 'metacom' and 'vegan' packages
  mat_Z <- aperm(Z, c(2,1))
  mat_Y <- aperm(Y, c(2,1))
  # Convert to zeros or ones for Y (instead of binomial output)
  mat_Y <- ifelse(mat_Y > 0, 1, 0)
  
  # Storage for metacommunity metrics:
  meta_Z <- NULL
  meta_Y <- NULL
  CoherZY <- array(0, dim=c(5, 2)) # 2 Columns: Z and Y
  TurnZY <- array(0, dim=c(5, 2)) # 2 Columns: Z and Y
  BoundZY <- array(0, dim=c(3, 2)) # 2 Coumns: Z and Y
  
  meta_Z <- tryCatch(Metacommunity(mat_Z, method="r1", sims=1000, allow.empty=T),
                     error = function(e){"ERROR"})
  meta_Y <- tryCatch(Metacommunity(mat_Y, method="r1", sims=1000, allow.empty=T),
                     error = function(e){"ERROR"})
  
  ####################################
  ### DETERMINE IF ERRORS OCCURRED ###
  ####################################
  
  ### If an error occurred, record this and return to the next iteration of "i"
  ### If an error has not occurred, complete the rest of the code, including the 
  ### Bayesian estimates of Z (posterior) and metacommunity structure.
  
  if(meta_Z == "ERROR" | meta_Y == "ERROR"){
    
    ERROR[i] <- TRUE
    Converge[i] <- "ERROR"
    
  }else { 
    ERROR[i] <- FALSE
    # Run the rest of the code:
    ############################################
    ######## BAYESIAN ESTIMATES FOR Z ##########
    ############################################
    # Update the data:
    X = X
    Y = Y
    K = ncol(Y)
    N = nrow(Y)
    J = J
    # current working directory:
    current_dir <- paste(getwd())
    
    # Make a new working directory (unique per node) name:
    wd.name <- paste("iter", runif(1,1,10), sep="_")
    system(paste("mkdir", wd.name, sep=" "))
    
    # Copy JAGS script files to new working directory:
    system(paste("cp", "./jags.script.cmd", paste("./", wd.name, "/jags.script.cmd", sep="")))
    system(paste("cp", "./OccMod_SingleYear.txt", paste("./", wd.name, "/OccMod_SingleYear.txt", sep="")))
    
    # Now change working directory to new one:
    setwd(paste("./", wd.name, sep=""))
    
    # Write a file of the data:
    dump(c("X", "Y", "K", "N", "J"), file="data.R")
    
    # Set initial parameters:
    # Z values (unobserved)
    zinit <- NULL
    zinit <- ifelse(Y > 0, 1, 0)
    z = zinit
    
    # Write a file of the inits:
    dump("z", file="inits.R")
    
    # Run the JAGS script file:
    system("jags jags.script.cmd")
    
    # read the coda files:
    ch1 <- NULL; ch2 <- NULL; ch3 <- NULL;
    ch1 <- read.coda("CODAchain1.txt", "CODAindex.txt", quiet=T)
    ch2 <- read.coda("CODAchain2.txt", "CODAindex.txt", quiet=T)
    ch3 <- read.coda("CODAchain3.txt", "CODAindex.txt", quiet=T)
    
    out <- NULL
    out <- mcmc.list(list(ch1, ch2, ch3))
    
    # Check for convergence for these two parameters:
    post.p <- NULL
    Rhat_p <- NULL
    
    post.p <- ggs(out, family="p")
    Rhat_p <- get_Rhat(post.p)
    
    if(any(Rhat_p > 1.1, na.rm=T)){
      Converge[i] <- FALSE }else{ Converge[i] <- TRUE}
    
    # Re-set working directory:
    setwd("../")
    
    # Delete all newly created files and directory:
    system(paste("rm -r", wd.name, sep=" "))
  }
}

