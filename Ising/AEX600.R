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


setwd('/Users/bokgyeongkang/work/Diagnostics/Ising')
#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
source("RFtns.R")

load("data/sim020Ising.RData")
# load("data/sim043Ising.RData")


#========================================================================
# for ACD
#========================================================================
# number of particles
dd = 600


aux.par = rep(0, dd)

### step 1. Conduct Liang's Fractional DMH to get auxiliary parameters
# N1 = 20000
# 
# sourceCpp("RcppFtns.cpp")
# ptm = proc.time()[3]
# FLiang = IsingFDMH(N1, 1, MPLE, 0.1, X)  # multiply 0.5 in c code
# timeFLiang = proc.time()[3] - ptm
# save(FLiang, timeFLiang, N1, file = 'ACD/AEX/simFDMH.RData')

load('ACD/AEX/simFDMH.RData')

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


### step 2. Run AEX
burnin = 1000

# Niter = 20*100000 + burnin
Niter = 20*50000 + burnin
Numaux = 10000*dd
t0 = 25000
neighbor = 20
sigma = 0.1
cycle = 1


sourceCpp("RcppFtns.cpp")

ptm = proc.time()[3]
res = IsingAEX(Niter, Numaux, cycle, t0, neighbor, aux.par, dist.aux.par, MPLE, sigma, X)
AEX = res$par
timeAEX = proc.time()[3] - ptm

save(burnin, Numaux, t0, neighbor, timeAEX, AEX, file = paste0('ACD/AEX/simAEX', dd, '.RData'))


