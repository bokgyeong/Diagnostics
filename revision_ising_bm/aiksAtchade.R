rm(list=ls())
library(Rcpp); library(RcppArmadillo)
library(doParallel); library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#================================================================
# Call data and functions
#================================================================
dd = 10
# dd = 20
# dd = 50
# dd = 100
# dd = 200
# dd = 400


load(paste0('ACDAIKS/Atchade/simAppxAtchade', dd, '.RData'))


#================================================================
# AIKS
#================================================================
nchains = length(appx)
score = foreach(i = 1:nchains, .combine = c) %do% { appx[[i]][,1] }
nrep = foreach(i = 1:nchains, .combine = c) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j] == Atchade[,i])) }
score = rep(score, nrep)
Atchade = foreach(i = 1:nchains, .combine = c) %do% { Atchade[,i] }

nn = 120000
# nn = 100000
# nn = 50000
indice = seq(20, nn, by = 20) # thinning
score = score[indice]
Atchade = Atchade[indice]
niter = length(Atchade)

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

# -----------------------------
nprocs = 19

ptm = proc.time()[3]
aiksAtchade = sum( apply(pIsingAIKS_woWeight(matrix(Atchade, ncol = 1), matrix(score, ncol = 1), c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
timeaiksAtchade = proc.time()[3] - ptm

save(aiksAtchade, timeaiksAtchade, file = paste0('ACDAIKS/Atchade/n', nn, '/simAIKSAtchade', dd, '.RData'))
# -----------------------------
