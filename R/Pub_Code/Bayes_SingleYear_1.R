# Single-year occupancy model and Bayesian analysis
# Author: JR Mihaljevic
# January 2014


############################################
######## BAYESIAN ESTIMATES FOR Z ##########
############################################

library(rjags)
library(ggmcmc)
setwd("~/Documents/Thesis Research/Occ.Model.EMS/R")

# Storage for output
post.z <- list()

for(i in 1:iter){
  
  jags_d <- list(X=X[, i],
                 Y=Yobs_mats[[i]],
                 K=ncol(Yobs_mats[[i]]),
                 N=nrow(Yobs_mats[[i]]),
                 J=J)
  
  # Set initial parameters:
  # Z values (unobserved)
  
  zinit <- ifelse(Yobs_mats[[1]] > 0, 1, 0)
  
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
  
  #Store output of z_ij values. 
  post.z[[i]] <- ggs(out, family="z")
  
}


####################################################
######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
####################################################

# Create a list with 'iter' elements.
# Each element is another list with the 1000 Z posterior estimates for each iter
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
