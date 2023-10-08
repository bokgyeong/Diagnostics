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

mm = c(1:5, 20)
n = 200000
# n = 250000
N = 200000

for(Nin in mm){
  load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))

  # ACD -----------------------------------------------------
  nchains = length(appx)

  D = foreach(i = 1:nchains) %do% { appx[[i]][,3:5] }
  nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Liang[[i]][,1])) }
  D = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:3, function(j) rep(D[[i]][,j], nrep[[i]])) }
  
  source('RFtns.R')
  ptm = proc.time()[3]
  acdLiang = ACD_bm(D[1:n,], N)
  timeacdLiang = proc.time()[3] - ptm

  cat('m =', Nin, ': ', acdLiang, '(', qchisq(0.99, 3), ')\n')

  save(acdLiang, timeacdLiang, file = paste0('ACDAIKS/Liang/n', n, '/simACDLiang', Nin, '.RData'))
}


