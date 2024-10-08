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
summax = 3
# summax = 5
# summax = 8
# summax = 9
# summax = 10
# summax = 20
# summax = 30
# summax = 100

appx = list()
start = 1; timeappx = 0
# load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))
# start = length(appx) + 1


load('data/bids2.RData')
load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', summax, '.RData'))
# load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', summax, 'v2.RData'))


#========================================================================
# biased but consistent approximation
#========================================================================
Trunc = Trunc[-(1:burn),]
niter = nrow(Trunc)

# nsets = 5
# nsets = 10
# nsets = 25
nsets = 20
# nsets = 35
Trunc = sapply(1:nsets, function(i) Trunc[((i-1) * (niter/nsets) + 1):(i * niter/nsets),], simplify = F)
th = sapply(1:nsets, function(i) unique(Trunc[[i]]), simplify = F)
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
  
  save(Trunc, th, nth, Sx, appx, timeappx, nu, file = paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))
  # save(Trunc, th, nth, Sx, appx, timeappx, nu, file = paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, 'v2.RData'))
}
