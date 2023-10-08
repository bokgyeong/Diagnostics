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

mm = c(100, 400, 450, 490, 500, 510)

# n = 200000
n = 120000
# n = 100000
# n = 50000

for(Nin in mm){
  load(paste('ACDAIKS/AEX/simAppxAEX', Nin, '.RData', sep = ''))

  # ACD -----------------------------------------------------
  nchains = length(appx)
  
  D = foreach(i = 1:nchains) %do% { appx[[i]][,2] }
  
  nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j]==AEX[, i])) }
  D = foreach(i = 1:nchains, .combine = c) %do% { rep(D[[i]], nrep[[i]]) }
  
  source('RFtns.R')
  ptm = proc.time()[3]
  acdAEX = ACD_bm(D[1:n])
  timeacdAEX = proc.time()[3] - ptm
  
  cat('d =', Nin, ':', acdAEX, '\n')
  
  save(acdAEX, timeacdAEX, file = paste0('ACDAIKS/AEX/n', n, '/simACDAEX', Nin, '.RData'))
}

qchisq(0.99, 1)
