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
# nprocs = 9


#========================================================================
# Load data
#========================================================================
summax = 3
# summax = 8
# summax = 9
# summax = 10
# summax = 30



load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))



# ------------------------------------------------------------------------------
# AIKS and p-value
# ------------------------------------------------------------------------------
p = ncol(Trunc[[1]])
nchains = length(appx)
score = foreach(i = 1:nchains) %do% { appx[[i]][,1:p] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1] == Trunc[[i]][,1])) }
score = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:p, function(j) rep(score[[i]][,j], nrep[[i]])) }
Trunc = foreach(i = 1:nchains, .combine = rbind) %do% { Trunc[[i]] }

# nn = 400000
nn = 300000
indice = seq(100, nn, by = 100) # thinning
score = score[indice,]
Trunc = Trunc[indice,]
niter = nrow(Trunc)


c = 1
beta = -1/2
nb = 1000

sourceCpp("RcppFtns.cpp")


# -----------------------------
ptm = proc.time()[3]
aiksTrunc = sum( apply(pAIKS_woWeight(Trunc, score, c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
timeaiksTrunc = proc.time()[3] - ptm
save(aiksTrunc, timeaiksTrunc, file = paste0('ACDAIKS/Trunc_knu/n', nn, '/bids2AIKSTrunc', summax, '.RData'))
# -----------------------------


# -----------------------------
# aiksTrunc = sum( apply(pAIKS_woWeight(Trunc, score, c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
# save(aiksTrunc, file = paste0('ACDAIKS/Trunc_knu/n', nn, '/bids2AIKSTrunc', summax, '.RData'))
# 
# 
# cl = parallel::makeCluster(nprocs, type = "PSOCK")
# doParallel::registerDoParallel(cl)
# 
# nbs = unique(c(0, seq(nprocs, nb, by = nprocs), nb)); aiksBootTrunc = c(); timeaiksTrunc = 0
# # load(paste0('ACDAIKS/Trunc_knu/n', nn, '/bids2AIKSTrunc', summax, '.RData'))
# # nbs = unique(c(seq(length(aiksBootTrunc), nb, by = nprocs), nb))
# 
# for(i in 2:length(nbs)){
#   ptm = proc.time()[3]
#   dummy = foreach(m = 1:(nbs[i]-nbs[i-1]), .combine = 'c', .packages = 'Rcpp', .noexport = c('AIKS', 'pAIKS')) %dopar% {
#     sourceCpp("RcppFtns.cpp")
#     qn = rmultinom(1, niter, rep(1/niter, niter))/niter - 1/niter
#     sum( apply(pAIKS(Trunc, score, qn, c, beta, niter, 1), 2, sum) )
#   }
#   timeaiksTrunc = timeaiksTrunc + proc.time()[3] - ptm
#   aiksBootTrunc = c(aiksBootTrunc, dummy)
#   save(aiksTrunc, aiksBootTrunc, timeaiksTrunc, file = paste0('ACDAIKS/Trunc_knu/n', nn, '/bids2AIKSTrunc', summax, '.RData'))
# }
# -----------------------------
