rm(list=ls())
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(snow)
library(doParallel)
library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#========================================================================
# Load data
#========================================================================
# appx = list()
# start = 1; timeappx = 0
load('ACDAIKS/Murray_knu/bids2AppxMurray.RData')
start = length(appx) + 1

load('data/bids2.RData')
load('ACDAIKS/Murray_knu/bids2Murray.RData')


#========================================================================
# biased but consistent approximation
#========================================================================
Murray = Murray[-(1:burn),]
niter = nrow(Murray)

# nsets = 5
nsets = 25
Murray = sapply(1:nsets, function(i) Murray[((i-1) * (niter/nsets) + 1):(i * niter/nsets),], simplify = F)
th = sapply(1:nsets, function(i) unique(Murray[[i]]), simplify = F)
nth = sapply(1:nsets, function(i) nrow(th[[i]]))

Sx = t(X) %*% y
nu = 1.75214

### simulation by Gibbs sampler
N = 200000 


nprocs = 19
# nprocs = 18
# nprocs = 7
mp_type = "PSOCK"
cl = parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)



for(j in start:nsets){
  
  ptm = proc.time()[3]
  appx[[j]] = foreach(i = 1:(nth[j]), .combine = 'rbind', .packages = "Rcpp") %dopar% {
    source('RFtns.R')
    sourceCpp("RcppFtns.cpp")
    
    theta = th[[j]][i,]
    Sy = rCOMP2_parallel_knu(X, theta, nu, N, 1)
    
    Uhat = f_Uhat_knu(Sx, Sy, nu)
    dhat = f_dhat_knu(Sx, Sy, nu, Uhat)
    
    c(Uhat, dhat)
  }
  timeappx = timeappx + proc.time()[3] - ptm
  
  save(Murray, th, nth, Sx, appx, timeappx, nu, file = 'ACDAIKS/Murray_knu/bids2AppxMurray.RData')
}

