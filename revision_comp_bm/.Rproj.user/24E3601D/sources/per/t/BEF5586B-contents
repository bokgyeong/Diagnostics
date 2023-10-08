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



#================================================================
# ACD using biased but consistent approximations
#================================================================
al = 0.01

# n = 400000
n = 300000
N = 200000
mm = c(3, 8, 9, 10, 30)

# ACD -----------------------------------------------------
for(summax in mm){
  
  load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))

  nsets = length(appx)
  p = ncol(Trunc[[1]])
  
  nrep = foreach(i = 1:nsets) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Trunc[[i]][,1] & th[[i]][j,p]==Trunc[[i]][,p])) }
  dhat = foreach(i = 1:nsets) %do% { appx[[i]][, (p + 1):(p + p*(p+1)/2)] }
  dhat = foreach(i = 1:nsets, .combine = rbind) %do% { sapply(1:(p*(p+1)/2), function(j) rep(dhat[[i]][,j], nrep[[i]])) }
  
  ptm = proc.time()[3]
  acdTrunc = ACD_bm(dhat[1:n,], N)
  timeacdTrunc = proc.time()[3] - ptm
  
  
  ### Look at each parameter
  cat('k =', summax, ':', acdTrunc, '(', qchisq(1-al, p*(p+1)/2), ')\n')

  save(acdTrunc, timeacdTrunc, file = paste0('ACDAIKS/Trunc_knu/n', n, '/bids2ACDTrunc', summax, '.RData'))
}

