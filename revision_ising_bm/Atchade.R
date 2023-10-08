rm(list=ls())
library(fields)
library(lattice)
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(xtable)
library(snow)
library(doParallel)
library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
load("data/sim020Ising.RData")


#========================================================================
# ALR
#========================================================================
# dd = 10
# dd = 20
# dd = 50
# dd = 100
# dd = 200
# dd = 400

aux.par = rep(0, dd)

### step 1. Conduct Liang's Fractional DMH to get auxiliary parameters
load('ACDAIKS/AEX/simFDMH.RData')

FLiang = FLiang[-(1:500)]                               # burn in 500
stand =  (FLiang-min(FLiang))/(max(FLiang)-min(FLiang)) # standardized
stand = unique(stand)                                   # only take unique components
dmat = rdist(stand)                                     # distance mat

# choose auxiliary parameters through min max procedure
ind = 1; A = 1; Ac = 2:length(stand)
aux.par[1] = stand[ind]

ind = which.max( dmat[,A] )
A = c(A,ind)
Ac = Ac[-which(Ac==ind)]
aux.par[2] = stand[ind]

for(i in 3:dd){
  dummy = max( apply( dmat[,A] , 1, min )[Ac] )
  ind = which(dmat[,A] == dummy, arr.ind = T)[1]
  A = c(A,ind)
  Ac = Ac[-which(Ac==ind)]
  aux.par[i] = stand[ind]
}

dist.aux.par = rdist(aux.par)  # distance matrix for aux.par (for standardized version)
aux.par = (max(FLiang)-min(FLiang))*aux.par + min(FLiang)


### step 2. Run ALR
burnin = 1000
Niter = 20*100000 + burnin

sourceCpp("RcppFtns.cpp")

ptm = proc.time()[3]
Atchade = IsingAtchade(Niter, 1, MPLE, 0.1, aux.par, X)
timeAtchade = proc.time()[3] - ptm
save(burnin, Atchade, timeAtchade, file = paste0('ACDAIKS/Atchade/simAtchade', dd, '.RData'))


