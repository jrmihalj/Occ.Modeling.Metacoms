joes_metacom_f <- function(p_mean=0.9, iter=1, z.iter=1000){
  
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
  # iter <- 1000
  
  # Storage (Outputs):
  Structure.ZY <- array("A", dim=c(iter, 2)) # 2 Columns: Z and Y
  Structure.Zpost <- array("A", dim=c(iter, z.iter)) # 1000 Z posteriors
  
  ##########################
  ###  START SIMULATION  ###
  ##########################
  
  for(i in 1:iter){
    
    # Establish basic parameters for the simulation:
    N <- 12 # Number of species
    K <- 75 # Number of sites
    J <- 4  # Number of sampling replicates per site
    
    # Fixed occurrence across species, in logit space
    b0 <- rep(Logit(0.6), N) # Change if desired
    # Detection probabilities:
    mu_p <- NULL
    mu_p <- Logit(p_mean)
    sd_p <- 0.75 #St. Dev. fixed for all simulations
    
    #### Draw covariate effects and values ####
    b.spp <- NULL; X <- NULL; lp0 <- NULL; p0 <- NULL
    
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
      
      Structure.ZY[i, 1:2] <- rep("ERROR", 2)
      Structure.Zpost[i, 1:z.iter] <- rep("ERROR", times=z.iter)
      
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
      # Update the data:
      X = X
      Y = Y
      K = ncol(Y)
      N = nrow(Y)
      J = J
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
      ch1 <- read.coda("CODAchain1.txt", "CODAindex.txt", quiet=T)
      ch2 <- read.coda("CODAchain2.txt", "CODAindex.txt", quiet=T)
      ch3 <- read.coda("CODAchain3.txt", "CODAindex.txt", quiet=T)
      
      out <- mcmc.list(list(ch1, ch2, ch3))
      
      # Check for convergence species-specific prob. detection (n=N):
      post.p <- NULL
      Rhat_p <- NULL
      
      post.p <- ggs(out, family="p")
      Rhat_p <- get_Rhat(post.p)
      
      # If not converged, update model until convergence is achieved
      if(any(Rhat_p > 1.1)){
        Structure.Zpost[i, 1:z.iter] <- rep("Not_Converged", times=z.iter)
      }else{
        #Store z_ij values. 
        post.z <- NULL
        post.z <- ggs(out, family="z")
        
        ####################################################
        ######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
        ####################################################
        
        for(j in 1:z.iter){ # 1000 draws from the posterior
          
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
          
          meta_Zpost <- tryCatch(Metacommunity(mat_Zpost, method="r1", sims=1000, allow.empty=T),
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
  # return:
  # Structure.ZY
  # Structure.Zpost
  results <- list(Structure.ZY = Structure.ZY,
                  Structure.Zpost = Structure.Zpost)
  return(results)
}