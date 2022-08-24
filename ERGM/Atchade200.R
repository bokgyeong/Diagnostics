rm(list=ls())
library(fields)
library(coda)
library(ergm)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(snow)
library(doParallel)
library(foreach)

setwd('/Users/bokgyeongkang/work/Diagnostics/ERGM')
#======================================================================
# call data and functions
#======================================================================
load("data/simERGM.RData")


# ======================================================================
# Atchade
# ======================================================================
dd = 200

aux.par = matrix(0, dd, 2)


### step 1. Conduct Liang's Fractional DMH to get particles
# N1 = 50000
# 
# sourceCpp("RcppFtns.cpp")
# 
# ptm = proc.time()[3]
# FLiang = ergmFDMH(X, diag(0.0025,2), matrix(hat,1), N1, 1)  # multiply 0.5 in c code
# timeFLiang = proc.time()[3] - ptm
# save(N1, FLiang, timeFLiang, file = 'ACD/Atchade/simFLiang.RData')

load('ACD/Atchade/simFLiang.RData')

FLiang = FLiang[-(1:500),]      # burn in 500
stand = matrix(0,nrow(FLiang),2)        # standardized
for(i in 1:2){ stand[,i] = (FLiang[,i]-min(FLiang[,i]))/(max(FLiang[,i])-min(FLiang[,i])) }
stand = unique(stand)
dmat = rdist(stand)             # distance matrix

# choose auxiliary parameters through min max procedure
ind = 1; A = 1; Ac = 2:dim(stand)[1]
aux.par[1,] = stand[ind,]

ind = which.max( dmat[,A] )
A = c(A,ind)
Ac = Ac[-which(Ac==ind)]
aux.par[2,] = stand[ind,]

for(i in 3:dd){
  dummy = max( apply( dmat[,A] , 1, min )[Ac] )
  ind = which(dmat[,A] == dummy, arr.ind = T)[1]
  A = c(A,ind)
  Ac = Ac[-which(Ac==ind)]
  aux.par[i,] = stand[ind,]
}

dist.aux.par = rdist(aux.par) # distance matrix for aux.par (for standardized version)
for(i in 1:2){ aux.par[,i] = (max(FLiang[,i])-min(FLiang[,i]))*aux.par[,i] + min(FLiang[,i]) }



### step 2. Run ALR
burnin = 1000
Niter = 25 * 75000 + burnin
COV = cov(aux.par)
cycle = 1


sourceCpp("RcppFtns.cpp")

ptm = proc.time()[3]
Atchade = ergmAtchade(Niter, cycle, COV, matrix(hat,1), aux.par, X)
timeAtchade = proc.time()[3] - ptm

save(burnin, X, stat, Atchade, timeAtchade, file = paste0('ACD/Atchade/simAtchade', dd, '.RData'))


