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
# Call functions
#========================================================================
source("RFtns.R")
sourceCpp("RcppFtns.cpp")



set.seed(1)
#========================================================================
# Simulate Data
#========================================================================
N = 30
# X = ProppWilson(N, N, 0.20)   # assume it is true data (you can change true parameter value)
X = ProppWilson(N, N, 0.43)
z = optim(0.1,pseudoLik,control=list(fnscale=-1), method = "BFGS")
MPLE = z$par

# save(X, N, MPLE, file = "simIsing.RData")
save(X, N, MPLE, file = "simIsing043.RData")


set.seed(1)
# X = ProppWilson(10, 10, 0.20)
X = ProppWilson(10, 10, 0.43)
Sy = pAuxStatPerf(10, 10, 0.20, 1000, 1)
ts.plot(Sy)
mean(Sy)
abline(h=Energy(X), col="red"); Energy(X)

