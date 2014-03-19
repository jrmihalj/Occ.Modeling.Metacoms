library(parallel)
library(Rmpi)
library(snow, lib.loc = "/home/majo3748/Rpackages")
library(permute, lib.loc = "/home/majo3748/Rpackages")
library(vegan, lib.loc = "/home/majo3748/Rpackages")
library(metacom, lib.loc = "/home/majo3748/Rpackages")

#Make the cluster
cl <- makeCluster( mpi.universe.size(), type="MPI" )

#Check the cluster - we should get a response from each worker
clusterCall( cl, function() Sys.info()[c("nodename","machine")])
#mpi.remote.exec(paste(Sys.info()[c("nodename")],"checking in as",mpi.comm.rank(),"of",mpi.comm.size()))


#Evaluate a function on each worker
joes_f <- function(iter=1){
  require(metacom)
  require(reshape2)
  Logit <- function(x){
    log(x) - log(1 - x)
  }

  AntiLogit <- function(x){
    exp(x) / (exp(x) + 1)
  }

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

  # Storage for Structure identities
  Structure <- array("A", dim=c(nrow(mat.eff), nrow(mat.val), iter))
  
  # simulate incidence matrices
  for (i in 1:nrow(mat.eff)){
    for(j in 1:nrow(mat.val)){
      for(it in 1:iter){      
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
      
        ### CHECK MATRIX FOR ROW/COLUMN SUMS ###
        # Check column sums
        if(any(colSums(Z) == 0) == TRUE){
          Z <- Z[, -which(colSums(Z) == 0)]
        }
        
        # Check row sums
        if(any(rowSums(Z) == 0) == TRUE){
          Z <- Z[-which(rowSums(Z) == 0), ]
        }

        # CALCULATE AND STORE METACOMMUNITY METRICS #        
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
          
          ###   DETERMINE METACOMMUNITY STRUCTURE   ####
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
  return(Structure.melted)
}

xx <- clusterCall(cl, joes_f, iter=10)
save(xx,file='~/test.rdata')

#Shutdown
stopCluster(cl)
mpi.quit()
