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
library(batchmeans)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("RFtns.R")

#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
load("data/sim020Ising.RData")


#========================================================================
# Liang
#========================================================================
# Nin = 1
# # Nin = 2
# # Nin = 3
# # Nin = 4
# # Nin = 5
# # Nin = 6
# # Nin = 7
# 
# burnin = 1000
# Nout = 20*100000 + burnin
# 
# sourceCpp("RcppFtns.cpp")
# 
# ptm = proc.time()[3]
# Liang = IsingDMH(Nout, Nin, MPLE, 0.1, X)
# timeLiang = proc.time()[3] - ptm
# 
# save(burnin, Liang, timeLiang, file = paste0('ACDAIKS_final/Liang/simLiang', Nin, '.RData'))




#========================================================================
# For gelman rubin diagnostic
#========================================================================
burnin = 1000
Nout = 120000 + burnin
Nin = rep(1:3, each = 3)
start = c(0.1, 0.3, 0.5, 0.1, 0.3, 0.5, 0.1, 0.3, 0.5)


### set up for paralellization
# nprocs = 19
# # mp_type = "MPI"
# mp_type = "PSOCK"
# cl = parallel::makeCluster(nprocs, type=mp_type)
# doParallel::registerDoParallel(cl)


### double MH
ptm = proc.time()[3]
# Liang = foreach(i = 1:length(Nin), .packages = "Rcpp", .noexport = c("Energe", "Gibb", "IsingDMH")) %dopar% {
Liang = foreach(i = 1:length(Nin), .packages = "Rcpp", .noexport = c("Energe", "Gibb", "IsingDMH")) %do% {
  sourceCpp("RcppFtns.cpp")
  IsingDMH(Nout, Nin[i], start[i], 0.1, X)
}
timeLiang = proc.time()[3] - ptm
save(burnin, Nin, Liang, start, timeLiang, file = paste0('mcmcDiag/simLiang.RData'))


