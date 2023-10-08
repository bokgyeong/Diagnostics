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
summax = 8
# summax = 9
# summax = 10
# summax = 20
# summax = 30
# summax = 100


load('data/bids2.RData')
Sx = t(X) %*% y
nu = 1.75214

### simulation by Gibbs sampler
N = 200000 
# nn = 400000
nn = 300000


# set up for parallelization ----------
cl = parallel::makeCluster(19, type = "PSOCK")
doParallel::registerDoParallel(cl)
# ------------------------------------

# outers = unique(c(0, seq(15, 100, by = 15), 100))

# start = 1
# start = 2
# start = 3
# start = 4
# start = 5
# start = 6
# start = 7

# for(repi in (outers[start]+1):(outers[start+1])){
for(repi in c(49, 56)){

  
  # appx = c(); timeappx = 0; starti = 1
  load(paste0('rep100/Trunc', summax, '/rep', repi, '_bids2AppxTrunc', summax, '.RData'))
  starti = nrow(appx) + 1
  
  load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', summax, '.RData'))
  Trunc = Trunc[burn+1:nn,]
  niter = nrow(Trunc)
  th = unique(Trunc)
  nth = nrow(th)
  
  
  ptm = proc.time()[3]
  dummy = foreach(i = starti:nth, .combine = 'rbind', .packages = "Rcpp") %dopar% {
    source('RFtns.R')
    sourceCpp("RcppFtns.cpp")
    
    theta = th[i,]
    Sy = rCOMP2_parallel_knu(X, theta, nu, N, 1)
    
    Uhat = f_Uhat_knu(Sx, Sy, nu)
    dhat = f_dhat_knu(Sx, Sy, nu, Uhat)
    
    c(Uhat, dhat)
  }
  timeappx = timeappx + proc.time()[3] - ptm
  
  appx = rbind(appx, dummy)
  
  save(Trunc, th, nth, Sx, appx, timeappx, nu, file = paste0('rep100/Trunc', summax, '/rep', repi, '_bids2AppxTrunc', summax, '.RData'))
}


