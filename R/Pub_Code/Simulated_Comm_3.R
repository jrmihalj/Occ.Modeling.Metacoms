# Determining bias in assigning metacommunity structure due to detection error
# Does an occupancy model do a better job at assigning 'known' structure of simulated communities?
# Author: JR Mihaljevic
# March 2014

##########################################
##########       TO DO:       ############
##########################################

# 1. CHANGE THE NUMBER OF ITERATIONS, IF DESIRED (HOW MANY UNIQUE COMMUNITIES TO SIMULATE)
#     - DEFAULT AT 1000

# 2. SET YOUR WORKING DIRECTORY TO WHERE THE JAGS OCCUPANCY MODEL IS LOCATED

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
require(rjags) # and all dependencies
require(ggmcmc) # and all dependencies

#######################
#  SET UP SIMULATION  #
#######################

# Set number of iterations (How many unique Z and Y matrices to create?)
iter <- 1000

# Establish basic parameters for the simulation:
N <- 12 # Number of species
K <- 75 # Number of sites
J <- 4  # Number of sampling replicates per site

# Base-line occurrence probability
# Fixed occurrence across species, in logit space
b0 <- rep(Logit(0.6), N) # Change if desired

# Fixed detection probability across species.
# The simulation will cycle through each of these:
p0s <- c(1.00, 0.75, 0.50)

# Storage (Outputs):
Structure.ZY <- array("A", dim=c(iter, 2)) # 2 Columns: Z and Y
Structure.Zpost <- array("A", dim=c(iter, 1000)) # 1000 Z posteriors

##########################
###  START SIMULATION  ###
##########################

for(p in 1:length(p0s)){
  
  #### RUN THE ENTIRE CODE FOR EACH OF THE SPECIFIED DETECTION PROBABILITIES ####
  p0 <- rep(p0s[p], N)
  
  
  #### RUN THE BULK OF THE CODE ####
  for(i in 1:iter){
    #### Draw covariate effects and values ####
    b.spp <- NULL
    b.spp <- rnorm(N, 0, runif(1, 0.5, 1)) # mean held at 0
    X <- NULL
    X <- rnorm(K, runif(1, 0, 2), runif(1, 2, 5))
    
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
    
    meta_Z <- tryCatch(Metacommunity(mat_Z, method="r1", sims=1000),
                       error = function(e){"ERROR"})
    meta_Y <- tryCatch(Metacommunity(mat_Y, method="r1", sims=1000),
                       error = function(e){"ERROR"})
    
    ####################################
    ### DETERMINE IF ERRORS OCCURRED ###
    ####################################
    
    ### If an error occurred, record this and return to the next iteration of "i"
    ### If an error has not occurred, complete the rest of the code, including the 
    ### Bayesian estimates of Z (posterior) and metacommunity structure.
    
    if(meta_Z == "ERROR" | meta_Y == "ERROR"){
      
      Structure.ZY[i, 1:2] <- rep("ERROR", 2)
      Structure.Zpost[i, 1:1000] <- rep("ERROR", 1000)
      
    }else { # Run the rest of the code:
      #############################################
      # CALCULATE AND STORE METACOMMUNITY METRICS #
      ############    FOR Z and Y    ##############
      #############################################
      
      CoherZY[, 1] <- as.numeric(as.character(meta_Z[[2]][1:5, ]))
      TurnZY[, 1] <- as.numeric(as.character(meta_Z[[3]][1:5, ]))
      
      CoherZY[, 2] <- as.numeric(as.character(meta_Y[[2]][1:5, ]))
      TurnZY[, 2] <- as.numeric(as.character(meta_Y[[3]][1:5, ]))
      
      for(g in 1:3){
        BoundZY[g, 1] <- meta_Z[[4]][1,g]
        BoundZY[g, 2] <- meta_Y[[4]][1,g]
      }
      
      #############################################
      ###   DETERMINE METACOMMUNITY STRUCTURE   ###
      ############    FOR Z and Y    ##############
      #############################################
      for(ZY in 1:2){
        
        if(CoherZY[3, ZY] > 0.05){
          Structure.ZY[i, ZY] <- "Random"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] > CoherZY[4, ZY]){
          Structure.ZY[i, ZY] <- "Checkerboard"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] > 0.05 & TurnZY[1, ZY] < TurnZY[4, ZY]){
          Structure.ZY[i, ZY] <- "Quasi-Nested"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] > 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[1, ZY] > 1 & BoundZY[2, ZY] <= 0.05){
          Structure.ZY[i, ZY] <- "Quasi-Clementsian"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] > 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[1, ZY] < 1 & BoundZY[2, ZY] <= 0.05){
          Structure.ZY[i, ZY] <- "Quasi-EvenSpaced"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] > 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[2, ZY] > 0.05){
          Structure.ZY[i, ZY] <- "Quasi-Gleasonian"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] <= 0.05 & TurnZY[1, ZY] < TurnZY[4, ZY]){
          Structure.ZY[i, ZY] <- "Nested"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] <= 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[1, ZY] > 1 & BoundZY[2, ZY] <= 0.05){
          Structure.ZY[i, ZY] <- "Clementsian"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] <= 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[1, ZY] < 1 & BoundZY[2, ZY] <= 0.05){
          Structure.ZY[i, ZY] <- "EvenSpaced"
        }
        
        if(CoherZY[3, ZY] <= 0.05 & CoherZY[1, ZY] < CoherZY[4, ZY] & 
             TurnZY[3, ZY] <= 0.05 & TurnZY[1, ZY] > TurnZY[4, ZY] &
             BoundZY[2, ZY] > 0.05){
          Structure.ZY[i, ZY] <- "Gleasonian"
        }
      }
      
      ############################################
      ######## BAYESIAN ESTIMATES FOR Z ##########
      ############################################
      jags_d <- NULL
      jags_d <- list(X=X,
                     Y=Y,
                     K=ncol(Y),
                     N=nrow(Y),
                     J=J)
      
      # Set initial parameters:
      # Z values (unobserved)
      zinit <- NULL
      zinit <- ifelse(Y > 0, 1, 0)
      
      # Start the model
      params <- c("lpsiMean", "lpsiSD", "z")
      
      mod <- NULL
      mod <- jags.model(file = "OccMod_SingleYear.txt", 
                        data = jags_d, n.chains = 3, n.adapt=1000,
                        inits = list(z=zinit))
      update(mod, n.iter=5000) # 5000 burn-in
      
      out <- NULL
      out <- coda.samples(mod, n.iter = 10000, variable.names = params, thin=10)
      
      # Check for convergence for these two parameters:
      post.psimean <- NULL
      post.psisd <- NULL
      Rhat_mean <- NULL
      Rhat_sd <- NULL
      
      post.psimean <- ggs(out, family="lpsiMean")
      post.psisd <- ggs(out, family="lpsiSD")
      Rhat_mean <- get_Rhat(post.psimean)
      Rhat_sd <- get_Rhat(post.psisd)
      
      # If not converged, update model until convergence is achieved
      if(Rhat_mean > 1.1 | Rhat_sd > 1.1){
        repeat{
          update(mod, n.iter=5000) # An extra 5,000 burn-in
          
          out <- NULL
          out <- coda.samples(mod, n.iter = 10000, variable.names = params, thin=10)
          
          post.psimean <- NULL
          post.psisd <- NULL
          Rhat_mean <- NULL
          Rhat_sd <- NULL
          post.psimean <- ggs(out, family="lpsiMean")
          post.psisd <- ggs(out, family="lpsiSD")
          Rhat_mean <- get_Rhat(post.psimean)
          Rhat_sd <- get_Rhat(post.psisd)
          
          if(Rhat_mean < 1.1 & Rhat_sd < 1.1) {break}
        }
      }
      
      #Store z_ij values. 
      post.z <- NULL
      post.z <- ggs(out, family="z")
      
      ####################################################
      ######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
      ####################################################
      
      for(j in 1:1000){ # 1000 draws from the posterior
        
        # Randomly choose from which chain the sample will originate:
        chain <- NULL
        chain <- sample(c(1:3), 1)
        
        # Create a subset vector of the values for all z[n, k]:
        subset <- NULL
        subset <- subset(post.z, Chain==chain & Iteration==paste(j))$value
        # Store this vector as a matrix
        Zpost <- NULL
        Zpost <- matrix(subset, nrow=nrow(Y), ncol=ncol(Y), byrow=F)
        
        
        # Check for zero sum columns/rows
        if(any(colSums(Zpost) == 0) == TRUE){
          Zpost <- Zpost[, -which(colSums(Zpost) == 0)]
        }
        if(any(rowSums(Zpost) == 0) == TRUE){
          Zpost <- Zpost[-which(rowSums(Zpost) == 0), ]
        }
        
        #############################################
        # CALCULATE AND STORE METACOMMUNITY METRICS #
        ##########    FOR EACH Z_POST    ############
        #############################################
        
        mat_Zpost <- NULL
        mat_Zpost <- aperm(Zpost, c(2,1)) # Transpose
        
        # Store all metacommunity metrics
        meta_Zpost <- NULL
        Coher.Zpost <- NULL
        Turn.Zpost <- NULL
        Bound.Zpost <- NULL
        
        meta_Zpost <- tryCatch(Metacommunity(mat_Zpost, method="r1", sims=1000),
                               error = function(e){"ERROR"})
        
        ####################################
        ### DETERMINE IF ERRORS OCCURRED ###
        ####################################
        
        if(meta_Zpost == "ERROR"){
          Structure.Zpost[i, j] <- "ERROR"
        }else {
          #############################################
          # CALCULATE AND STORE METACOMMUNITY METRICS #
          ############    FOR Z.post     ##############
          #############################################
          
          Coher.Zpost <- as.numeric(as.character(meta_Zpost[[2]][1:5, ]))
          Turn.Zpost <- as.numeric(as.character(meta_Zpost[[3]][1:5, ]))
          for(q in 1:3){
            Bound.Zpost[q] <- meta[[4]][1,q]
          }
          
          #############################################
          ###   DETERMINE METACOMMUNITY STRUCTURE   ###
          ############    FOR Z.post     ##############        
          #############################################
          
          if(Coher.Zpost[3] > 0.05){
            Structure.Zpost[i, j] <- "Random"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] > Coher.Zpost[4]){
            Structure.Zpost[i, j] <- "Checkerboard"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] > 0.05 & Turn.Zpost[1] < Turn.Zpost[4]){
            Structure.Zpost[i, j] <- "Quasi-Nested"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] > 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[1] > 1 & Bound.Zpost[2] <= 0.05){
            Structure.Zpost[i, j] <- "Quasi-Clementsian"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] > 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[1] < 1 & Bound.Zpost[2] <= 0.05){
            Structure.Zpost[i, j] <- "Quasi-EvenSpaced"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] > 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[2] > 0.05){
            Structure.Zpost[i, j] <- "Quasi-Gleasonian"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] <= 0.05 & Turn.Zpost[1] < Turn.Zpost[4]){
            Structure.Zpost[i, j] <- "Nested"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] <= 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[1] > 1 & Bound.Zpost[2] <= 0.05){
            Structure.Zpost[i, j] <- "Clementsian"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] <= 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[1] < 1 & Bound.Zpost[2] <= 0.05){
            Structure.Zpost[i, j] <- "EvenSpaced"
          }
          
          if(Coher.Zpost[3] <= 0.05 & Coher.Zpost[1] < Coher.Zpost[4] & 
               Turn.Zpost[3] <= 0.05 & Turn.Zpost[1] > Turn.Zpost[4] &
               Bound.Zpost[2] > 0.05){
            Structure.Zpost[i, j] <- "Gleasonian"
          }
          
        }
        
      }
      
    }
    
  }
  
}


####################
#  END SIMULATION  #
####################

# Outputs:

# Structure.ZY
# Structure.Zpost

