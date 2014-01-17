# Function to plot an ordinated incidence matrix as a heat map in ggplot2
# August 2013
# Author: JR Mihaljevic

#Input an ORDINATED matrix, and the appropriate x and y labels
#Outputs a plot of the matrix in black and white. 
Matrix_Plot <- function(matrix, xlab="Parasite Species", ylab="Host Individuals"){
  require(reshape2)
  require(ggplot2)
  
  ord_mat <- matrix
  #Need to number sites so the graph doesn't rearrange matrix
  rownames(ord_mat) <- as.factor(c(1:nrow(ord_mat)))
  Names <- NULL
  Names <- colnames(ord_mat)
  colnames(ord_mat) <- as.factor(1:length(Names))
  
  # melt() changes the data into a data.frame that ggplot2 can read for the heat map
  dataframe <- melt(ord_mat)
  colnames(dataframe) <- c("Site", "Species", "Pres")
  
  # Plot matrix in ggplot2
  incidence <- ggplot(dataframe, aes(x=Species, y=Site))+
    geom_tile(aes(fill=factor(Pres)), color=element_blank())+
    labs(x=xlab, y=ylab)+
    scale_fill_discrete(h=c(0, 180), h.start=0, l=c(100,0), c=0)+
    theme_classic()+
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    #scale_x_discrete(limits=c(1:length(Names)), labels=Names)+ #Optional to label species
    #theme(axis.text.x=element_text(angle=60, vjust=0.5))+
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(legend.position="none")+
    scale_y_reverse()
  
  L <- list(plot = incidence)
}

# Example:
library(metacom)
data(testmatrices)
test_ord <- OrderMatrix(testmatrices[[7]])

Fig4a <- Matrix_Plot(test_ord, xlab="Species", ylab="Sites")

quartz(height=6, width=3)
print(Fig4a$plot)

