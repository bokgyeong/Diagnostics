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
# ACD 
#================================================================

# mm = 1:7
mm = 1:5

# n = 200000
n = 120000
# n = 100000
# n = 50000

for(Nin in mm){
  load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))

  # ACD -----------------------------------------------------
  nchains = length(appx)

  U = foreach(i = 1:nchains) %do% { appx[[i]][,1] }
  D = foreach(i = 1:nchains) %do% { appx[[i]][,2] }

  nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j]==Liang[, i])) }
  U = foreach(i = 1:nchains, .combine = c) %do% { rep(U[[i]], nrep[[i]]) }
  D = foreach(i = 1:nchains, .combine = c) %do% { rep(D[[i]], nrep[[i]]) }

  source('RFtns.R')
  ptm = proc.time()[3]
  acdLiang = ACD_bm(D[1:n])
  timeacdLiang = proc.time()[3] - ptm

  cat('m =', Nin, ':', acdLiang, '\n')

  save(acdLiang, timeacdLiang, file = paste0('ACDAIKS/Liang/n', n, '/simACDLiang', Nin, '.RData'))
}

qchisq(0.99, 1)


plot(th[[1]][order(th[[1]])], D[[1]][order(th[[1]])])

