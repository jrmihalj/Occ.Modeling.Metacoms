library(parallel)
library(Rmpi)
library(snow, lib.loc = "/home/majo3748/Rpackages")
library(permute, lib.loc = "/home/majo3748/Rpackages")
library(vegan, lib.loc = "/home/majo3748/Rpackages")
library(metacom, lib.loc = "/home/majo3748/Rpackages")
library(reshape2)

#Make the cluster
cl <- makeCluster( mpi.universe.size(), type="MPI" )

#Check the cluster - we should get a response from each worker
clusterCall( cl, function() Sys.info()[c("nodename","machine")])
#mpi.remote.exec(paste(Sys.info()[c("nodename")],"checking in as",mpi.comm.rank(),"of",mpi.comm.size()))

#Evaluate a function on each worker
source("joes_f.r")

xx <- clusterCall(cl, joes_f, iter=1)
save(xx,file='test.rdata')

#Shutdown
stopCluster(cl)
mpi.quit()
