library(doMC)
registerDoMC(4)

joes_f <- function(iter=1){
  # Permuations of covariate effects:
  # Normal (0, 0.5)
  # Normal (0, 1)
  # Uniform (-1, 1)
  # Uniform (-3, 3)
  
  # Permutations of covariate values (for each covariate effect):
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
  b0 <- rep(Logit(0.6), N) # Fixed occurrence across species, in logit space
  
  # Establish species-specific detection probabilities
  p0 <- rep(0.7, N) # Fixed detection prob.
  
  
  #############################################
  ### SIMULATE COMMUNITY INCIDENCE MATRICES ###
  #############################################
  
  # Set number of iterations for each set of parameter values
  total.iter <- iter*nrow(mat.eff)*nrow(mat.val) # Total number of simulations
  
  # Storage
  b.spp <- vector() # species-specific covariate effects
  X <- vector() # site-specific covariate values
  
  z <- array(0, dim = c(N, K)) # 'True' occupancy states
  Y <- array(0, dim = c(N, K, nrow(mat.eff), nrow(mat.val), iter)) # This will hold simulated observations
  lpsi <- array(0, dim = c(N, K)) # This will hold the occupancy probabilities
  psi <- array(0, dim = c(N, K))
  
  
  for (i in 1:nrow(mat.eff)){
    for(j in 1:nrow(mat.val)){
      for(it in 1:iter){
        
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
  
  #############################################
  ####  STORE SEPARATE MATRICES IN A LIST  ####
  #############################################
  
  All_mats <- list()
  All_mats_names <- paste("Mat_", 1:total.iter, sep="")
  count <- 1
  for (it in 1:iter){
    for(i in 1:nrow(mat.eff)){
      for(j in 1:nrow(mat.val)){
        All_mats[[count]] <- Y[, , i, j, it]
        names(All_mats)[count] <- All_mats_names[i]
        count <- count + 1
      }
    }
  }
  
  #############################################
  ### CHECK EACH MATRIX FOR ROW/COLUMN SUMS ###
  #############################################
  
  # First, check row and column sums.
  # Second, transpose so that rows are sites and columns are species
  # This is necessary for the 'metacom' package, which uses 'vegan' functions. 
  # 'vegan' requires a site-by-species matrix. 
  
  for(i in 1:total.iter){
    # Check column sums
    if(any(colSums(All_mats[[i]]) == 0) == TRUE){
      print(paste(which(colSums(All_mats[[i]]) == 0),"This site of Mat", i, "sums to zero", sep=" "))
      All_mats[[i]] <- All_mats[[i]][, -which(colSums(All_mats[[i]]) == 0)]
    }
    # Check row sums
    if(any(rowSums(All_mats[[i]]) == 0) == TRUE){
      print(paste(which(rowSums(All_mats[[i]]) == 0), ":","This species of Mat", i, "sums to zero", sep=" "))
      All_mats[[i]] <- All_mats[[i]][-which(rowSums(All_mats[[i]]) == 0), ]
    }
    # Transpose
    All_mats[[i]] <- aperm(All_mats[[i]], c(2, 1))
    # Convert to pres/abs
    All_mats[[i]] <- ifelse(All_mats[[i]] > 0, 1, 0)
  }
  
  #############################################
  # CALCULATE AND STORE METACOMMUNITY METRICS #
  ##########    FOR EACH MATRIX    ############
  #############################################
  
  library(metacom)
  
  #Storage of metacommunity metrics
  Coher <- array(0, dim=c(total.iter, 5))
  colnames(Coher) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
  Turn <- array(0, dim=c(total.iter, 5))
  colnames(Turn) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
  Bound <- array(0, dim=c(total.iter, 3))
  colnames(Bound) <- c("index", "pval", "df")
  
  for(i in 1:total.iter){
    meta <- NULL
    meta <- Metacommunity(All_mats[[i]], method="r1", sims=1000)
    Coher[i, ] <- as.numeric(as.character(meta[[2]][1:5, ]))
    Turn[i, ] <- as.numeric(as.character(meta[[3]][1:5, ]))
    for(j in 1:3){
      Bound[i, j] <- meta[[4]][1,j]
    }
  }
  
  ##############################################
  ###   DETERMINE METACOMMUNITY STRUCTURE   ####
  ###########    FOR EACH MATRIX    ############
  ##############################################
  
  # Store output:
  Structure <- vector()
  
  for(i in 1:total.iter){
    if(Coher[i, 3] > 0.05){
      Structure[i] <- "Random"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] > Coher[i, 4]){
      Structure[i] <- "Checkerboard"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] > 0.05 & Turn[i, 1] < Turn[i, 4]){
      Structure[i] <- "Quasi-Nested"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 1] > 1 & Bound[i, 2] <= 0.05){
      Structure[i] <- "Quasi-Clementsian"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 1] < 1 & Bound[i, 2] <= 0.05){
      Structure[i] <- "Quasi-EvenSpaced"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 2] > 0.05){
      Structure[i] <- "Quasi-Gleasonian"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] <= 0.05 & Turn[i, 1] < Turn[i, 4]){
      Structure[i] <- "Nested"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 1] > 1 & Bound[i, 2] <= 0.05){
      Structure[i] <- "Clementsian"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 1] < 1 & Bound[i, 2] <= 0.05){
      Structure[i] <- "EvenSpaced"
    }
    
    if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
         Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
         Bound[i, 2] > 0.05){
      Structure[i] <- "Gleasonian"
    }
  }
  
  
  #########################################
  ###   STORE ALL INPUTS AND RESULTS   ####
  #########################################
  
  Cov.Eff.Dist <- rep(rep(paste("Dist", 1:nrow(mat.eff), sep=""), each=nrow(mat.val)), times=iter)
  Cov.Val.Dist <- rep(paste("Dist", 1:nrow(mat.val), sep=""), times=nrow(mat.eff)*iter)
  
  Output <- data.frame(Cov.Eff.Dist, Cov.Val.Dist, Structure)
  return(Output)
}

out <- foreach(i=1:getDoParWorkers(), .combine=rbind) %dopar% {
  subout <- joes_f(iter=125)
  subout
}
write.csv(out, file="sim_dist_results.csv")
