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



###############################
# Max ran the simulations on JANUS
# Output .RData file with object 'xx'
# 'xx': list of 84 elements, each a dataframe

Z_data <- NULL
for(i in 1:length(xx)){
  Z_data <- rbind(Z_data, xx[[i]])
}
colnames(Z_data) <- c("Ef.dist", "Val.dist","Iter","Structure")

# Get rid of "ERROR" ones
remove <- which(Z_data$Structure=="ERROR") #139
Z_data <- Z_data[-remove, ]
Z_data$Structure <- factor(Z_data$Structure)
levels(Z_data$Structure)

# Make a table to summarize results:
# I want the count of each structure observed for each combination of 
# Cov.Eff.Dist and Cov.Val.Dist
library(reshape2)
casted <- dcast(Z_data, Ef.dist + Val.dist ~ Structure)
Z_melt <- melt(casted, id.vars=c("Ef.dist", "Val.dist"))
colnames(Z_melt) <- c("Ef.dist", "Val.dist", "Structure", "Count")
Z_melt$Ef.dist <- as.factor(Z_melt$Ef.dist)
Z_melt$Val.dist <- as.factor(Z_melt$Val.dist)

#Change Ef.dist labels so I can order them on the graph:
levels(Z_melt$Ef.dist)[levels(Z_melt$Ef.dist)=="1"] <- "N(0,0.5)"
levels(Z_melt$Ef.dist)[levels(Z_melt$Ef.dist)=="2"] <- "N(0,1)"
levels(Z_melt$Ef.dist)[levels(Z_melt$Ef.dist)=="3"] <- "U(-3,0)"

# Graph the frequencies of each structure:
# First, set up a color-blind-friendly color palette:
cbPalette <- c("#999999", "#E69F00", "#CC79A7", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#56B4E9", "#000000")

library(ggplot2)
frac <- ggplot(Z_melt, aes(x=Val.dist, y=Count, fill=Structure))+
  geom_bar(stat="identity")+
  theme_classic()+
  facet_grid(Ef.dist~.)+
  scale_y_continuous(breaks=c(0, 500, 1000))+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values=cbPalette, breaks=c("Random", "Clementsian", "Quasi-Clementsian",
                                               "Gleasonian", "Quasi-Gleasonian", "Nested", 
                                               "Quasi-Nested", "EvenSpaced", 
                                               "Checkerboard"))+
  scale_x_discrete(limits=c("1","5","2","6","3","7","4","8"),
                   labels=c("N(0,1)","U(-1,1)","N(0,5)","U(-4,4)","N(2,1)","U(0,3)","N(5,2)","U(0,10)"))+
  theme(strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(angle=90, hjust=.5))

quartz(height=10, width=10)
print(frac) 
