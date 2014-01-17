# Construct the Z-mat posterior and calculate/visualize EMS posteriors
# Author: JR Mihaljevic
# January 2014

####################################################
######## CONSTRUCT Z FROM BAYESIAN OUTPUT ##########
####################################################

iter <- 100 # number of samples to draw from the posterior 

# Each chain has 300 observations (n.iter=3000, thinned by 10)
# Choose which iterations will be used from the posterior
samples <- sample(c(1:300), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(nrow(Y), ncol(Y), iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.out, Chain==chain & Iteration==samples[i])$value
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

for(i in 1:iter){
  meta <- NULL
  meta <- Metacommunity(z.post[, , i], method="r1", sims=1000, allow.empty=T)
  z.ord[, , i] <- meta[[1]]
  Coher[i, ] <- as.numeric(as.character(meta[[2]][1:5, ]))
  Turn[i, ] <- as.numeric(as.character(meta[[3]][1:5, ]))
  for(j in 1:3){
    Bound[i, j] <- meta[[4]][1,j]
  }
}

# Generate heat map of ordinated matrices:
quartz(height=6, width=3)
print(Matrix_HeatMap(z.ord, ylab="Sites"))

# Generate density plots of posterior metacommunity metrics:
Coher.plot <- ggplot(data.frame(Coher), aes(x=z))+
  geom_histogram()+
  labs(x="Coherence z-score", y="histogram")+
  theme_classic()

Turn.plot <- ggplot(data.frame(Turn), aes(x=z))+
  geom_histogram()+
  labs(x="Turnover z-score", y="")+
  theme_classic()

Bound.plot <- ggplot(data.frame(Bound), aes(x=index))+
  geom_histogram()+
  labs(x="Morista's I", y="")+
  theme_classic()


library(gridExtra)
grob <- arrangeGrob(Coher.plot, Turn.plot, Bound.plot, nrow=1)
quartz(height=5, width=9)
print(grob)



