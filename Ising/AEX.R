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
source("RFtns.R")

#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
load("data/sim020Ising.RData")
# load("data/sim043Ising.RData")


#========================================================================
# AEX
#========================================================================
# nchains = 100
# 
# # number of particles
# dd = 30
# # dd = 60
# # dd = 100
# 
# aux.par = rep(0, dd)
# 
# ### step 1. Conduct Liang's Fractional DMH to get auxiliary parameters
# # N1 = 20000
# # 
# # sourceCpp("RcppFtns.cpp")
# # set.seed(1016)
# # ptm = proc.time()[3]
# # FLiang = IsingFDMH(N1, 10, MPLE, 0.1, X)  # multiply 0.5 in c code
# # timeFLiang = proc.time()[3] - ptm
# # save(FLiang, timeFLiang, N1, file = 'AEX/simFDMH.RData')
# 
# 
# load('AEX/simFDMH.RData')
# 
# FLiang = FLiang[-(1:1000)]                               # burn in 1000
# stand =  (FLiang-min(FLiang))/(max(FLiang)-min(FLiang)) # standardized
# stand = unique(stand)                                   # only take unique components
# dmat = rdist(stand)                                     # distance mat
# 
# # choose auxiliary parameters through min max procedure
# ind = 1; A = 1; Ac = 2:length(stand)
# aux.par[1] = stand[ind]
# 
# ind = which.max( dmat[,A] )
# A = c(A,ind)
# Ac = Ac[-which(Ac==ind)]
# aux.par[2] = stand[ind]
# 
# for(i in 3:dd){
#   dummy = max( apply( dmat[,A] , 1, min )[Ac] )
#   ind = which(dmat[,A] == dummy, arr.ind = T)[1]
#   A = c(A,ind)
#   Ac = Ac[-which(Ac==ind)]
#   aux.par[i] = stand[ind]
# }
# 
# dist.aux.par = rdist(aux.par)  # distance matrix for aux.par (for standardized version)
# aux.par = (max(FLiang)-min(FLiang))*aux.par + min(FLiang)
# 
# 
# ### step 2. Run AEX
# # burnin = 1000000
# burnin = 100000
# 
# Niter = 1000*20 + burnin
# Numaux = 10000*dd
# t0 = 25000
# neighbor = 20
# sigma = 0.1
# 
# cycle = 1
# # cycle = 5
# # cycle = 10
# 
# 
# nprocs = 19
# mp_type = "PSOCK"
# cl = parallel::makeCluster(nprocs, type=mp_type)
# doParallel::registerDoParallel(cl)
# 
# 
# ptm = proc.time()[3]
# AEX = foreach(i = 1:nchains, .combine = 'cbind', .packages = "Rcpp") %dopar% {
#   source("RFtns.R")
#   sourceCpp("RcppFtns.cpp")
#   
#   # res = IsingAEX(Niter, Numaux, cycle, t0, neighbor, aux.par, dist.aux.par, MPLE,sigma, X)
#   # res$par[-(1:burnin)]
#   
#   IsingAEX2(Niter, Numaux, cycle, t0, neighbor, aux.par, dist.aux.par, MPLE,sigma, X)
# }
# timeAEX = proc.time()[3] - ptm
# 
# # save(burnin, Numaux, t0, neighbor, timeAEX, AEX, file = paste0('AEX/simAEX', dd, 'm', cycle, '.RData'))
# save(burnin, Numaux, t0, neighbor, timeAEX, AEX, file = paste0('AEX/simAEX2d', dd, 'm', cycle, '.RData'))



#========================================================================
# AEX2
#========================================================================
# nchains = 5
# 
# # number of particles
# dd = 50
# # dd = 100
# # dd = 200
# # dd = 400
# 
# aux.par = rep(0, dd)
# 
# ### step 1. Conduct Liang's Fractional DMH to get auxiliary parameters
# # N1 = 20000
# # 
# # sourceCpp("RcppFtns.cpp")
# # set.seed(1016)
# # ptm = proc.time()[3]
# # FLiang = IsingFDMH(N1, 10, MPLE, 0.1, X)  # multiply 0.5 in c code
# # timeFLiang = proc.time()[3] - ptm
# # save(FLiang, timeFLiang, N1, file = 'AEX/simFDMH.RData')
# 
# 
# load('AEX/simFDMH.RData')
# 
# FLiang = FLiang[-(1:500)]                               # burn in 500
# stand =  (FLiang-min(FLiang))/(max(FLiang)-min(FLiang)) # standardized
# stand = unique(stand)                                   # only take unique components
# dmat = rdist(stand)                                     # distance mat
# 
# # choose auxiliary parameters through min max procedure
# ind = 1; A = 1; Ac = 2:length(stand)
# aux.par[1] = stand[ind]
# 
# ind = which.max( dmat[,A] )
# A = c(A,ind)
# Ac = Ac[-which(Ac==ind)]
# aux.par[2] = stand[ind]
# 
# for(i in 3:dd){
#   dummy = max( apply( dmat[,A] , 1, min )[Ac] )
#   ind = which(dmat[,A] == dummy, arr.ind = T)[1]
#   A = c(A,ind)
#   Ac = Ac[-which(Ac==ind)]
#   aux.par[i] = stand[ind]
# }
# 
# dist.aux.par = rdist(aux.par)  # distance matrix for aux.par (for standardized version)
# aux.par = (max(FLiang)-min(FLiang))*aux.par + min(FLiang)
# 
# 
# ### step 2. Run AEX
# burnin = 1000
# 
# Niter = 20000*20 + burnin
# Numaux = 10000*dd
# t0 = 25000
# neighbor = 20
# sigma = 0.1
# cycle = 1
# 
# 
# nprocs = nchains
# mp_type = "PSOCK"
# cl = parallel::makeCluster(nprocs, type=mp_type)
# doParallel::registerDoParallel(cl)
# 
# 
# ptm = proc.time()[3]
# AEX = foreach(i = 1:nchains, .combine = 'cbind', .packages = "Rcpp") %dopar% {
#   source("RFtns.R")
#   sourceCpp("RcppFtns.cpp")
#   res = IsingAEX(Niter, Numaux, cycle, t0, neighbor, aux.par, dist.aux.par, MPLE,sigma, X)
#   res$par
# }
# timeAEX = proc.time()[3] - ptm
# 
# save(burnin, Numaux, t0, neighbor, timeAEX, AEX, file = paste0('AEX2/simAEX', dd, '.RData'))



#========================================================================
# for ACD
#========================================================================
# number of particles
# dd = 200
# dd = 300
# dd = 400
dd = 450
# dd = 500
# dd = 600


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


