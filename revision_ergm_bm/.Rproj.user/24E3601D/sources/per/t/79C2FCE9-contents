rm(list=ls())
library(coda)
library(ergm)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(snow)
library(doParallel)
library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#======================================================================
# call data and functions
#======================================================================
load("data/simERGM.RData")
# load("data/sim2ERGM.RData")


#======================================================================
# Liang
#======================================================================
# Nin = 1
Nin = 2
# Nin = 3
# Nin = 4
# Nin = 5
# Nin = 6
# Nin = 7
# Nin = 8
# Nin = 9
# Nin = 10
# Nin = 20

burnin = 1000
Niter = 25 * 100000 + burnin
th = matrix(hat, 1, 2)


sourceCpp("RcppFtns.cpp")
ptm = proc.time()[3]
Liang = ergmDMH(X, COV, th, Niter, Nin)[-1,]
timeLiang = proc.time()[3] - ptm

# save(burnin, X, stat, Liang, timeLiang, file = paste0('ACDAIKS/Liang/simLiang', Nin, '.RData'))




# ptm = proc.time()[3]
# Liang = ergmDMH(X, COV, th, 100, 3)[-1,]
# timeLiang = proc.time()[3] - ptm
# timeLiang/100*(75000*25+burnin)/60/60
