# Examples for  ESA talk 2014
setwd("~/Documents/Thesis Research/Occ.Model.EMS/R/ESA2014")

# Hypothetical metacom for title slide
hypomet <- read.csv(file="HypotheticalMetacom.csv", header=T)
hypomet <- as.matrix(hypomet)
colnames(hypomet) <- c(1:ncol(hypomet))
rownames(hypomet) <- c(1:nrow(hypomet))

quartz(height=5, width=6)
print(Matrix_Plot(hypomet, xlab="Species", ylab="Sites"))


# Reduced version:

smaller <- hypomet[1:65, 1:26]
quartz(height=5, width=6)
print(Matrix_Plot(smaller, xlab="Species", ylab="Sites"))

# Random small:
small_unord <- matrix(rbinom(65*26, 1, 0.3), nrow=65, ncol=26)
colnames(small_unord) <- c(1:ncol(small_unord))
rownames(small_unord) <- c(1:nrow(small_unord))

quartz(height=5, width=6)
print(Matrix_Plot(small_unord, xlab="Species", ylab="Sites"))

# Poke holes in the dataset
library(metacom)
small_holes <- smaller
for(i in 1:nrow(smaller)){
  for(j in 1:ncol(smaller)){
    if(smaller[i,j]==1){
      small_holes[i,j] <- rbinom(1,1,0.5)
    } else(small_holes[i,j]<-0)
  }
}

if(any(rowSums(small_holes) == 0) == TRUE){
  small_holes <- small_holes[-which(rowSums(small_holes) == 0), ]
}

Metacommunity(small_holes)
quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(small_holes), xlab="Species", ylab="Sites"))


# Gleasonian:
gleason <- read.csv(file="Gleasonian.csv", header=T)
gleason <- as.matrix(gleason)
colnames(gleason) <- 1:ncol(gleason)
rownames(gleason) <- 1:nrow(gleason)

quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(gleason), xlab="Species", ylab="Sites"))

gleason_holes <- gleason
for(i in 1:nrow(gleason)){
  for(j in 1:ncol(gleason)){
    if(gleason[i,j]==1){
      gleason_holes[i,j] <- rbinom(1,1,0.9)
    } else(gleason_holes[i,j]<-0)
  }
}

quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(gleason_holes), xlab="Species", ylab="Sites"))

# Nested:
nested <- read.csv(file="Nested.csv", header=T)
nested <- as.matrix(nested)
colnames(nested) <- 1:ncol(nested)
rownames(nested) <- 1:nrow(nested)

quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(nested), xlab="Species", ylab="Sites"))

nested_holes <- nested
for(i in 1:nrow(nested)){
  for(j in 1:ncol(nested)){
    if(nested[i,j]==1){
      nested_holes[i,j] <- rbinom(1,1,0.9)
    } else(nested_holes[i,j]<-0)
  }
}

if(any(rowSums(nested_holes) == 0) == TRUE){
  nested_holes <- nested_holes[-which(rowSums(nested_holes) == 0), ]
}
quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(nested_holes), xlab="Species", ylab="Sites"))


# Gleasonian with more detection error:
gleason_ImpDet <- gleason_holes
for(i in 1:nrow(gleason)){
  for(j in 1:ncol(gleason)){
    if(gleason[i,j]==1){
      gleason_ImpDet[i,j] <- rbinom(1,1,0.6)
    } else(gleason_ImpDet[i,j]<-0)
  }
}

quartz(height=5, width=6)
print(Matrix_Plot(OrderMatrix(gleason_ImpDet), xlab="Species", ylab="Sites"))

