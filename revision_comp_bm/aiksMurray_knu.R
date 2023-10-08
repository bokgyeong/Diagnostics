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


nprocs = 19


#========================================================================
# Load data
#========================================================================
load('ACDAIKS/Murray_knu/bids2AppxMurray.RData')


# ------------------------------------------------------------------------------
# AIKS and p-value
# ------------------------------------------------------------------------------
p = ncol(Murray[[1]])
nchains = length(appx)
score = foreach(i = 1:nchains) %do% { appx[[i]][,1:p] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1] == Murray[[i]][,1])) }
score = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:p, function(j) rep(score[[i]][,j], nrep[[i]])) }
Murray = foreach(i = 1:nchains, .combine = rbind) %do% { Murray[[i]] }

nn = 300000
# nn = 400000
indice = seq(100, nn, by = 100) # thinning
score = score[indice,]
Murray = Murray[indice,]
niter = nrow(Murray)

c = 1
beta = -1/2
nb = 1000


# -----------------------------
sourceCpp("RcppFtns.cpp")
ptm = proc.time()[3]
aiksMurray = sum( apply(pAIKS_woWeight(Murray, score, c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
timeaiksMurray = proc.time()[3] - ptm
save(aiksMurray, timeaiksMurray, file = paste0('ACDAIKS/Murray_knu/n', nn, '/bids2AIKSMurray.RData'))

# -----------------------------

