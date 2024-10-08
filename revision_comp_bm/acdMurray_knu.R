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
source('RFtns.R')


#================================================================
# ACD using biased but consistent approximations
#================================================================


n = 300000
# n = 400000
N = 200000

# ACD -----------------------------------------------------

load('ACDAIKS/Murray_knu/bids2AppxMurray.RData')

nsets = length(appx)
p = ncol(Murray[[1]])

nrep = foreach(i = 1:nsets) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Murray[[i]][,1] & th[[i]][j,p]==Murray[[i]][,p])) }
dhat = foreach(i = 1:nsets) %do% { appx[[i]][, (p + 1):(p + p*(p+1)/2)] }
dhat = foreach(i = 1:nsets, .combine = rbind) %do% { sapply(1:(p*(p+1)/2), function(j) rep(dhat[[i]][,j], nrep[[i]])) }

ptm = proc.time()[3]
acdMurray = ACD_bm(dhat[1:n,], N)
timeacdMurray = proc.time()[3] - ptm


### Look at each parameter
cat(acdMurray, '(', qchisq(0.99, p*(p+1)/2), ')\n')

save(acdMurray, timeacdMurray, file = paste0('ACDAIKS/Murray_knu/n', n, '/bids2ACDMurray.RData'))



