#!/bin/bash
#PBS -N mpi_test
#PBS -l walltime=15:00:00
#PBS -l nodes=21:ppn=12
#PBS -j oe

#We just asked for 1 node with 12 processors and 4 hours of wall time.
# (in the janus-short queue)

#The following will also load dependent modules, incl openmpi
module load R/r-3.0.1
module load userApps/jags

#Change directory to the working directory
cd $PBS_O_WORKDIR

orterun -np 1 -hostfile $PBS_NODEFILE R --vanilla -f mpi_f.r
