# JR Mihaljevic
# May 2014

# Analyzing/visualizing output of data, testing whether the occupancy model
# reduces bias in assigning metacommunity structure:

# Manipulating N (24 and 36)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Output (TEST)
# Using p_mean=0.50 and J=3

# 'xx' is the name of the Janus output file
str(xx)
# Structure of file is:
# List of 60 (cores)
# Each list has subset list of 2 objects:
#   - Structure.ZY: dim=c(2,2), 
#   - Structure.Zpost: dim=c(2, 500)

p.5_J3_N36 <- xx

#1: Concatinate the lists to form two objects
ZY_p.5_J3_N36 <- NULL # Z and Y matrices' metacommunity structures
Z.post_p.5_J3_N36 <- NULL # Z.posterior matrices' metacommunity structures

for(i in 1:120){
  ZY_p.5_J3_N36 <- rbind(ZY_p.5_J3_N36, p.5_J3_N36[[i]]$Structure.ZY)
  Z.post_p.5_J3_N36 <- rbind(Z.post_p.5_J3_N36, p.5_J3_N36[[i]]$Structure.Zpost)
}
#------------------------------------------------------------------------------#

#2: Determine if any had errors

unique(ZY_p.5_J3_N36[,2]) #Check Z and Y (no errors)
unique(Z.post_p.5_J3_N36[,1]) #None
#------------------------------------------------------------------------------#

#3: Determine the proportion of Y that match Z
YmatchZ <- length(which(ZY_p.5_J3_N36[,1]==ZY_p.5_J3_N36[,2])) / (120) #120 iterations
# 0.66

#------------------------------------------------------------------------------#

#4: Determine the proportion of Z.post that match Z for each iteration
ZPmatchZ <- NULL

for(i in 1:nrow(Z.post_p.5_J3_N36)){
  ZPmatchZ[i] <- length(which(Z.post_p.5_J3_N36[i, ]==ZY_p.5_J3_N36[i,1])) / 500
}

library(modeest)
mlv(ZPmatchZ, method = "mfv")
# Mode=1.00
mean(ZPmatchZ) #0.79
median(ZPmatchZ) #0.88

# Bootstrap estimated 95%CI of the median:
library(boot)
med <- function(y, indices) median(y[indices])
b <- boot(ZPmatchZ, med, 10000)
plot(b)
boot.ci(b, type="bca")
# ( 0.780,  0.895 ) N24
# ( 0.822,  0.931 ) N36
boot.ci(b, type="norm")
# ( 0.7996,  0.9172 ) N24
# ( 0.8264,  0.9271 ) N36

#------------------------------------------------------------------------------#

#5. Store as data.frame
Treat <- rep("p.5_J3_N36", 120)
df_Zpost_p.5_J3_N36 <- data.frame(Treat, ZPmatchZ)
#------------------------------------------------------------------------------#

#6. Plot the data (histogram)

library(ggplot2)

plot_p.5_J3_N36 <- ggplot(df_Zpost_p.5_J3_N36, aes(x=ZPmatchZ))+
  geom_histogram(binwidth=0.03, linetype=0)+
  theme_classic()+
#   labs(x=expression(paste("Prop. ", Z[post], " matching Z")),
#        y="Frequency")+
  labs(x="",y="")+
  scale_x_continuous(breaks=c(0, 0.5, 1))+
  scale_y_continuous(breaks=c(0, 15, 30))+
  theme(axis.text.y=element_text(angle=90, hjust=.5))+
  geom_vline(xintercept=YmatchZ, color="black")+
  geom_vline(xintercept=median(ZPmatchZ), color="gray")+
  geom_vline(xintercept=c(.82,.93), color="gray", linetype=2)

quartz(height=3, width=4.5)
print(plot_p.5_J3_N36)

#------------------------------------------------------------------------------#

# Plot the % correct by structure type:
df2_Zpost_p.5_J3_N36 <- data.frame(ZY_p.5_J3_N36[,1], ZPmatchZ)
colnames(df2_Zpost_p.5_J3_N36) <- c("Z", "ZPmatchZ")

plot2_p.5_J3_N36 <- ggplot(df2_Zpost_p.5_J3_N36, aes(x=Z, y=ZPmatchZ))+
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
print(plot2_p.5_J3_N36)
  
  
  
  