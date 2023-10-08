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


#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
load("data/sim020Ising.RData")
sourceCpp("RcppFtns.cpp")


#========================================================================
# Conduct Murray's Exchagne Algorithm
#========================================================================
burnin = 1000
Niter = 20*100000 + burnin


cpoint = seq(0, Niter, by = 1000)

# Murray = c(); timeMurray = 0; initial = MPLE; resume = 2
load('ACDAIKS_final/Murray/simMurray.RData')
resume = which(nrow(Murray) == cpoint) + 1

for(i in resume:length(cpoint)){
  
  ptm = proc.time()[3]
  res = IsingExchange(cpoint[i] - cpoint[i-1], initial, 0.1, X)
  timeMurray = proc.time()[3] - ptm + timeMurray
  
  Murray = rbind(Murray, res)
  initial = Murray[nrow(Murray),]
  save(X, Murray, timeMurray, initial, file = 'ACDAIKS_final/Murray/simMurray.RData')
}


