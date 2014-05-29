# JR Mihaljevic
# May 2014

# Analyzing/visualizing output of data, testing whether the occupancy model
# reduces bias in assigning metacommunity structure:

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Output (TEST)
# Using p_mean=0.50 and J=4 (will reduce to J=3 later)

# 'xx' is the name of the Janus output file
str(xx)
# Structure of file is:
# List of 252 (cores)
# Each list has subset list of 2 objects:
#   - Structure.ZY: dim=c(4,2), 
#   - Structure.Zpost: dim=c(4, 500)

p.5_J3 <- xx

#1: Concatinate the lists to form two objects
ZY_p.5_J3 <- NULL # Z and Y matrices' metacommunity structures
Z.post_p.5_J3 <- NULL # Z.posterior matrices' metacommunity structures

for(i in 1:252){
  ZY_p.5_J3 <- rbind(ZY_p.5_J3, p.5_J3[[i]]$Structure.ZY)
  Z.post_p.5_J3 <- rbind(Z.post_p.5_J3, p.5_J3[[i]]$Structure.Zpost)
}
#------------------------------------------------------------------------------#

#2: Determine if any had errors

unique(ZY_p.5_J3[,2]) #Check Z and Y (no errors)
unique(Z.post_p.5_J3[,1]) #None
#------------------------------------------------------------------------------#

#3: Determine the proportion of Y that match Z
YmatchZ <- length(which(ZY_p.5_J3[,1]==ZY_p.5_J3[,2])) / (252*4) #4 iterations/core
# 0.775

#------------------------------------------------------------------------------#

#4: Determine the proportion of Z.post that match Z for each iteration
ZPmatchZ <- NULL

for(i in 1:nrow(Z.post_p.5_J3)){
  ZPmatchZ[i] <- length(which(Z.post_p.5_J3[i, ]==ZY_p.5_J3[i,1])) / 500
}

library(modeest)
mlv(ZPmatchZ, method = "mfv")
# Mode=1.00
mean(ZPmatchZ) #0.78
median(ZPmatchZ) #0.934

# Bootstrap estimated 95%CI of the median:
library(boot)
med <- function(y, indices) median(y[indices])
b <- boot(ZPmatchZ, med, 10000)
plot(b)
boot.ci(b, type="bca")
# ( 0.913,  0.950 )
boot.ci(b, type="norm")
# ( 0.9161,  0.9537 )   

#------------------------------------------------------------------------------#

#5. Store as data.frame
Treat <- rep("p.5_J3", 4*252)
df_Zpost_p.5_J3 <- data.frame(Treat, ZPmatchZ)
#------------------------------------------------------------------------------#

#6. Plot the data (histogram)

library(ggplot2)

plot_p.5_J3 <- ggplot(df_Zpost_p.5_J3, aes(x=ZPmatchZ))+
  geom_histogram(binwidth=0.03, linetype=0)+
  theme_classic()+
#   labs(x=expression(paste("Prop. ", Z[post], " matching Z")),
#        y="Frequency")+
  labs(x="",y="")+
  scale_x_continuous(breaks=c(0, 0.5, 1))+
  scale_y_continuous(breaks=c(0, 75, 150))+
  theme(axis.text.y=element_text(angle=90, hjust=.5))+
  geom_vline(xintercept=YmatchZ, color="black")+
  geom_vline(xintercept=median(ZPmatchZ), color="gray")+
  geom_vline(xintercept=c(.784,.834), color="gray", linetype=2)

quartz(height=3, width=4.5)
print(plot_p.5_J3)

#------------------------------------------------------------------------------#

# Plot the % correct by structure type:
df2_Zpost_p.5_J3 <- data.frame(ZY_p.5_J3[,1], ZPmatchZ)
colnames(df2_Zpost_p.5_J3) <- c("Z", "ZPmatchZ")

plot2_p.5_J3 <- ggplot(df2_Zpost_p.5_J3, aes(x=Z, y=ZPmatchZ))+
  geom_point(size=.5, shape=1, alpha=0.9, color="lightgray")+
  stat_summary(fun.data="mean_cl_boot", conf.int=0.95,
               color="black", size=0.8)+
  labs(x="", y=expression(paste("Prop. ", Z[post], " matching Z")))+
  scale_y_continuous(breaks=c(0,.5,1), labels=c("0","0.5","1"))+
  theme_classic()+
#   theme(axis.text.y=element_text(size=10, angle=90, vjust=.5),
#         axis.text.x=element_text(angle, hjust=0.5))+
  coord_flip()

quartz(height=5, width=8)
print(plot2_p.5_J3)
  
  
  
  