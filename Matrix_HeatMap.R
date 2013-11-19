Matrix_HeatMap <- function(array, xlab="Species", ylab="Site"){
  require(reshape2)
  require(ggplot2)
  
  ord_mat <- array
  sum_mat <- array(0, dim=c(nrow(ord_mat), ncol(ord_mat)))
  for(i in 1:nrow(ord_mat)){
    for(j in 1:ncol(ord_mat)){
      sum_mat[i,j] <- sum(ord_mat[i, j, ])
    }
  }
  #Need to number sites so the graph doesn't rearrange matrix
  rownames(sum_mat) <- as.factor(c(1:nrow(sum_mat)))
  colnames(sum_mat) <- as.factor(c(1:ncol(sum_mat)))
  
  # melt() changes the data into a data.frame that ggplot2 can read for the heat map
  dataframe <- melt(sum_mat)
  colnames(dataframe) <- c("Site", "Species", "Pres")
  
  # Plot matrix in ggplot2
  incidence <- ggplot(dataframe, aes(x=Species, y=Site, fill=Pres))+
    geom_tile(color=element_blank())+
    labs(x=xlab, y=ylab)+
    scale_fill_gradient(low="yellow", high="red")+
    theme_classic()+
    theme(axis.text=element_blank(), axis.ticks=element_blank())+
    theme(legend.position="none")+
    scale_y_reverse()
  
  L <- list(plot = incidence)
}

quartz(height=6, width=3)
print(Matrix_HeatMap(Y.ord, xlab="Species", ylab="Sites"))