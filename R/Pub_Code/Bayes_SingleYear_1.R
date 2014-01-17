# Single-year occupancy model and Bayesian analysis
# Author: JR Mihaljevic
# January 2014

# Covariate effects for occurrence, but not for detection
# Model stored in document: OccMod_SingleYear.txt

# model {
#   #-------------------------------------------------------------------------------
#   #######################################
#   ############ ASSIGN PRIORS ############
#   #######################################
#   #-------------------------------------------------------------------------------
#   ## Metacommunity-wide average of parameters: psi, p:
#   
#   psiMean ~ dbeta(1,1) # Occurrence probability
#   pMean ~ dbeta(1,1) # Detection probability
#   
#   #-------------------------------------------------------------------------------
#   ## Manual logit transformation of mean parameter values:
#   
#   lpsiMean <- log(psiMean) - log(1-psiMean)
#   lpMean <- log(pMean) - log(1-pMean) 
#   
#   #-------------------------------------------------------------------------------
#   ## Precision estimates for each metacommunity average:
#   
#   lpsiSD ~ dunif(0,10)
#   lpsiPrec <- pow(lpsiSD,-2)
#   
#   lpSD ~ dunif(0,10)
#   lpPrec <- pow(lpSD,-2)
#   
#   #-------------------------------------------------------------------------------
#   ## Estimates for covariate effect means and precisions:
#   
#   bMean ~ dnorm(0,.001) # b = covariates for occurrence probability
#   bSD ~ dunif(0,10)
#   bPrec <- pow(bSD,-2)
#   
#   #-------------------------------------------------------------------------------
#   #######################################
#   ########## Likelihood Model ###########
#   #######################################
#   #-------------------------------------------------------------------------------
#   for(i in 1:N){
#     ## Species-specific baseline of occurrence prob.
#     b0[i] ~ dnorm(lpsiMean, lpsiPrec)T(-12,12) 
#     
#     ## Species-specific covariate effects on occurrence:
#     b[i] ~ dnorm(bMean, bPrec)
#     
#     ## Species-specific detection probability
#     lp[i] ~ dnorm(lpMean, lpPrec)T(-12,12)
#     p[i] <- 1/(1+exp(-lp[i])) # Manual anti-logit
#     
#     for(k in 1:K){
#       lpsi[i, k] <- b0[i] + b[i] * X[k]
#       psi[i, k] <- 1/(1+exp(-lpsi[i, k]))
#       
#       z[i, k] ~ dbern(psi[i, k])
#       Y[i, k] ~ dbinom(p[i] * z[i, k], J)
#     }  
#   }
#   # Close Model  
# }


####################################################
########## BAYESIAN OCCUPANCY MODELING #############
####################################################

jags_d <- list(X=X,
               Y=Y,
               K=ncol(Y),
               N=nrow(Y),
               J=J)

# Set initial parameters:
# Z values (unobserved)
zinit <- Y
zinit <- ifelse(zinit > 0, 1, 0)

# Start the model
library(rjags)

params <- c("z") # This will give me species-level estimates for each of these

mod <- jags.model(file = "OccMod_SingleYear.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))

out <- coda.samples(mod, n.iter = 3000, variable.names = params, thin=10)

# Store output:
library(ggmcmc)
post.out <- ggs(out)