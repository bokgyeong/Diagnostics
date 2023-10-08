rm(list=ls())
library(Rcpp); library(RcppArmadillo)
library(doParallel); library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#================================================================
# Call data and functions
#================================================================
dd = 100
# dd = 400
# dd = 450
# dd = 490
# dd = 500
# dd = 510

load(paste0('ACDAIKS/AEX/simAppxAEX', dd, '.RData'))



#================================================================
# AIKS
#================================================================
nchains = length(appx)
score = foreach(i = 1:nchains, .combine = c) %do% { appx[[i]][,1] }
nrep = foreach(i = 1:nchains, .combine = c) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j] == AEX[,i])) }
score = rep(score, nrep)
AEX = foreach(i = 1:nchains, .combine = c) %do% { AEX[,i] }

nn = 120000
# nn = 100000
# nn = 50000
indice = seq(20, nn, by = 20) # thinning
score = score[indice]
AEX = AEX[indice]
niter = length(AEX)

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

# -----------------------------
nprocs = 19

ptm = proc.time()[3]
aiksAEX = sum( apply(pIsingAIKS_woWeight(matrix(AEX, ncol = 1), matrix(score, ncol = 1), c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
timeaiksAEX = proc.time()[3] - ptm

save(aiksAEX, timeaiksAEX, file = paste0('ACDAIKS/AEX/n', nn, '/simAIKSAEX', dd, '.RData'))
# -----------------------------






