# Function for creating ordinated heat map for posterior Z-matrices
# Author: JR Mihaljevic
# January 2014

Matrix_HeatMap_NMDS <- function(matrix, xlab="Species", ylab="Site"){
  require(reshape2)
  require(ggplot2)
  
  ord_mat <- matrix
  
  #Need to number sites so the graph doesn't rearrange matrix
  rownames(ord_mat) <- as.factor(c(1:nrow(ord_mat)))
  Names <- NULL
  Names <- colnames(ord_mat)
  colnames(ord_mat) <- as.factor(c(1:ncol(ord_mat)))

  # Turn zeros into NA values (to make it white):
  ord_mat[which(ord_mat==0)] <- NA
  
  # melt() changes the data into a data.frame that ggplot2 can read for the heat map
  dataframe <- melt(ord_mat)
  colnames(dataframe) <- c("Site", "Species", "Pres")
  
  # Plot matrix in ggplot2
  incidence <- ggplot(dataframe, aes(x=Species, y=Site, fill=Pres))+
    geom_tile(color=element_blank())+
    labs(x=xlab, y=ylab)+
    scale_fill_gradient(low="white", high="black", na.value="transparent")+
    theme_classic()+
    theme(axis.text=element_blank(), axis.ticks=element_blank())+
#     theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     scale_x_discrete(limits=c(1:length(Names)), labels=Names)+
#     theme(axis.text.x=element_text(angle=60, vjust=0.5))+
    theme(legend.position="none")#+
    #scale_y_reverse()
  
  L <- list(plot = incidence)
}
