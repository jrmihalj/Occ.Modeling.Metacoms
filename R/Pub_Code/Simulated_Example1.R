# Simulated metacommunity example (first example presented in the text)
# Demonstrates how to integrate occupancy modeling with EMS methods
# Data stored in: Simulated_Example1.RData

# Author: JR Mihaljevic
# Date: March 2014

#####################################################################################
#-----------------------------------------------------------------------------------#
#####################################################################################

####################################################
##########  ESTABLISH USEFUL FUNCTIONS  ############
####################################################
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

##################################################
##########  LOAD REQUIRED R PACKAGES  ############
##################################################
require(metacom)
require(rjags)
require(ggmcmc)


####################################################
########## SIMULATION OF OCCUPANCY DATA ############
####################################################
# Set number of iterations (How many unique Z and Y matrices to create?)
# iter <- 1000

# Establish basic parameters for the simulation:
N <- 12 # Number of species
K <- 75 # Number of sites
J <- 3  # Number of sampling replicates per site

# Base-line occurrence probability
# Fixed occurrence across species, in logit space
b0 <- rep(Logit(0.6), N) # Change if desired

# Community-level detection probability (hyper-parameters)
p_mean <- .50
mu_p <- Logit(p_mean) # Logit space
sd_p <- 0.75 # standard dev.


######################
#  START SIMULATION  #
######################


b.spp <- NULL
b.spp <- rnorm(N, 0, 1)
X <- NULL
X <- rnorm(K, 0, 2)

lp0 <- rnorm(N, mu_p, sd_p) # Use if detection assigned by hyperparameters
p0 <- AntiLogit(lp0)

#### Occupancy states ####
# Storage:
Z <- array(0, dim = c(N, K)) # 'True' occupancy states
Y <- array(0, dim = c(N, K)) # Simulated observations
lpsi <- array(0, dim = c(N, K)) # occupancy probabilities (see below)
psi <- array(0, dim = c(N, K))
for(n in 1:N){
  for(k in 1:K){
    lpsi[n, k] <- b0[n] + b.spp[n] * X[k] # Covariate effects on occurrence 
    psi[n, k] <- AntiLogit(lpsi[n, k])
    
    Z[n, k] <- rbinom(1, 1, psi[n, k]) # True Occupancy
    Y[n, k] <- rbinom(1, J, p0[n] * Z[n, k]) # Observed Occupancy
  }
}



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
# Convert to zeros or ones (instead of binomial output)
#mat_Z <- ifelse(mat_Z > 0, 1, 0)
mat_Y <- ifelse(mat_Y > 0, 1, 0)

# Plot the ordinated Y and Z matrices:
quartz(height=6, width=3)
print(Matrix_Plot(OrderMatrix(mat_Z), xlab="", ylab=""))

rotate <- function(x) t(apply(x, 2, rev))
quartz(height=6, width=3)
print(Matrix_Plot(rotate(rotate(OrderMatrix(mat_Y))), xlab="", ylab=""))    

# Store all metacommunity metrics
meta_Z <- NULL
meta_Y <- NULL
meta_Z <- Metacommunity(mat_Z, method="r1", sims=1000)
meta_Y <- Metacommunity(mat_Y, method="r1", sims=1000)

# Plot the un-ordinated Y and Z matrices:
quartz(height=6, width=3)
print(Matrix_Plot(mat_Z, xlab="Species", ylab="Site"))

quartz(height=6, width=3)
print(Matrix_Plot(mat_Y, xlab="Species", ylab="Site"))

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
Rhat_mean <- get_Rhat(post.psimean) # RHAT MUST BE < 1.1
Rhat_sd <- get_Rhat(post.psisd)

#Store z_ij values. 
post.z <- NULL
post.z <- ggs(out, family="z")



######################################
#### SAMPLE FROM THE Z_POSTERIOR  ####
######################################

iter <- 1000 # number of samples to draw from the posterior 

# Each chain has 1000 observations (n.iter=10000, thinned by 10)
# Choose which iterations will be used from the posterior
samples <- sample(c(1:1000), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(nrow(Y), ncol(Y), iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.z, Chain==chain & Iteration==samples[i])$value
  # Store this vector as a matrix
  z.mat <- NULL
  z.mat <- matrix(subset, nrow=nrow(Y), ncol=ncol(Y), byrow=F)
  # Add matrix to the array
  z.post[, , i] <- z.mat
}

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
######## CALCULATE METACOMMUNITY METRICS  ##########
####################################################

# Make sites rows and species columns:
z.post <- aperm(z.post, c(2, 1, 3))

# Ordinate all the z matrices:

library(metacom)
# Store the Ordinated Matrices:
z.ord <- array(0, dim=c(ncol(Y), nrow(Y), iter))
# Store the statistics for Coherence and Turnover.
Coher <- array(0, dim=c(iter, 5))
colnames(Coher) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
Turn <- array(0, dim=c(iter, 5))
colnames(Turn) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
Bound <- array(0, dim=c(iter, 3))
colnames(Bound) <- c("index", "pval", "df")

pb <- txtProgressBar(min = 0, max = iter, style = 3)
for(i in 1:iter){
  meta <- NULL
  meta <- Metacommunity(z.post[, , i], method="r1", sims=1000)
  z.ord[, , i] <- meta[[1]]
  Coher[i, ] <- as.numeric(as.character(meta[[2]][1:5, ]))
  Turn[i, ] <- as.numeric(as.character(meta[[3]][1:5, ]))
  for(j in 1:3){
    Bound[i, j] <- meta[[4]][1,j]
  }
  setTxtProgressBar(pb, i)
}
# Determine Structure for each:
Structure <- NULL

for(i in 1:iter){
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

###################################
#### POSTERIORS OF EMS METRICS ####
###################################
library(ggplot2)

Coher.plot <- ggplot(data.frame(Coher), aes(x=z))+
  geom_histogram(binwidth=0.3)+
  geom_histogram(binwidth=0.3, data=data.frame(Coher[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  geom_vline(xintercept=1.96, linetype=5)+ # Significance cut-off
  scale_y_continuous(limits=c(0, 120), breaks=c(0, 60, 120))

quartz(height=3.5, width=3.5)
print(Coher.plot)
  
Turn.plot <- ggplot(data.frame(Turn), aes(x=z))+
  geom_histogram(binwidth=0.3)+
  geom_histogram(binwidth=0.3, data=data.frame(Turn[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  geom_vline(xintercept=c(-1.96), linetype=5)+
  scale_x_continuous(breaks=c(-9, -3, 3))+
  scale_y_continuous(limits=c(0, 90), breaks=c(0, 45, 90))

quartz(height=3.5, width=3.5)
print(Turn.plot)

Bound.plot <- ggplot(data.frame(Bound), aes(x=index))+
  geom_histogram(binwidth=0.02)+
  geom_histogram(binwidth=0.02, data=data.frame(Bound[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  scale_x_continuous(breaks=c(1, 1.3, 1.6))+
  scale_y_continuous(limits=c(0, 90), breaks=c(0, 45, 90))

quartz(height=3.5, width=3.5)
print(Bound.plot)

####################################
#### DISTRIBUTION OF STRUCTURES ####
####################################

head(Structure)
unique(Structure)

Structure2 <- Structure

# Change the order to something more logical/ideal for plotting
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Clementsian"), "b_Quasi-Clem")
Structure2 <- replace(Structure2, which(Structure2=="Clementsian"), "a_Clem")
Structure2 <- replace(Structure2, which(Structure2=="Gleasonian"), "c_Gleas")
Structure2 <- replace(Structure2, which(Structure2=="Random"), "f_Rand")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Nested"), "e_Quasi-Nest")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Gleasonian"), "d_Quasi_Gleas")

Str.Counts <- as.data.frame(table(Structure2))
Str.Counts$ID <- rep("A", nrow(Str.Counts)) # Need an ID var for the bar plot

Str.Counts
 

frac <- ggplot(Str.Counts, aes(x=ID, y=Freq, fill=Structure2))+
  geom_bar(stat="identity")+
  theme_classic()+
  labs(x="", y="")+
  scale_fill_grey(breaks=paste(Str.Counts$Structure),
                    labels=c("Clementsian: 65.3%", "Quasi-Clementsian: 6.7%", "Gleasonian: 3.9%",
                             "Quasi-Gleasonian: 0.1%", "Quasi-Nested: 0.3%", "Random: 23.7%"))+
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(angle=90, hjust=0.5), 
        axis.text.x=element_blank(), legend.title=element_blank())

quartz(height=5, width=4)
print(frac)

#####################################
####   HEAT MAP OF Z-POSTERIOR   ####
#####################################

Z_heat <- mat.or.vec(nr=K, nc=N)
for(i in 1:1000){
  Z_heat <- Z_heat + z.post[,,i]
}
Z_heat <- Z_heat/1000 # Now each cell is a probability of occurrence averaged over z.post
rownames(Z_heat) <- as.factor(1:K)
colnames(Z_heat) <- as.factor(1:N)
# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_heat.CCA <- decorana(Z_heat, ira=0) # CHECK CONVERGENCE, ETC.
Z_heat.CCA.sites <- Z_heat.CCA$rproj[, 1]
Z_heat.CCA.spp <- Z_heat.CCA$cproj[, 1]

Z_heat_ord <- Z_heat[order(Z_heat.CCA.sites, decreasing = FALSE), 
                     order(Z_heat.CCA.spp, decreasing = FALSE)]

quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_heat_ord, xlab="", ylab=""))
#################################################################################
#### LOOK AT PROB DETECTION AND SPECIES-SPECIFIC COVARIATE EFFECT POSTERIORS ####
#################################################################################
library(rjags)
library(ggmcmc)

mod2 <- jags.model(file = "OccMod_SingleYear.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))
update(mod2, n.iter=5000)
out2 <- coda.samples(mod2, n.iter = 10000, variable.names = c("p", "b"), thin=10)

post.p <- NULL
post.p <- ggs(out2, family="p")
Rhat_p <- get_Rhat(post.p)
head(post.p)

post.b <- NULL
post.b <- ggs(out2, family="b")
Rhat_b <- get_Rhat(post.b)
head(post.b)

###########################      
#######  DET. PROB  ####### 
###########################  
library(mcmcplots)

quartz(height=5, width=7)
caterplot(out2, parms="p", col="black", val.lim=c(0, 1))
caterpoints(p0)


##############################      
#######  COV. EFFECTS  ####### 
############################## 

quartz(height=5, width=7)
caterplot(out2, parms="b", col="black")
caterpoints(b.spp)
