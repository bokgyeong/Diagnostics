rm(list=ls())
library(Rcpp); library(RcppArmadillo)
library(doParallel); library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


nprocs = 19
# nprocs = 11



#================================================================
# Call data and functions
#================================================================
# dd = 1
# dd = 2
# dd = 3
# dd = 4
# dd = 5
dd = 20

load(paste0('ACDAIKS/Liang/simAppxLiang', dd, '.RData'))


# ------------------------------------------------------------------------------
# AIKS and p-value
# ------------------------------------------------------------------------------
p = ncol(Liang[[1]])
nchains = length(appx)
score = foreach(i = 1:nchains) %do% { appx[[i]][,1:p] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1] == Liang[[i]][,1])) }
score = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:p, function(j) rep(score[[i]][,j], nrep[[i]])) }
Liang = foreach(i = 1:nchains, .combine = rbind) %do% { Liang[[i]] }

# nn = 200000
nn = 250000
indice = seq(25, nn, by = 25) # thinning
score = score[indice,]
Liang = Liang[indice,]
niter = nrow(Liang)

c = 1
beta = -1/2
nb = 1000
# nb = 2000

sourceCpp("RcppFtns.cpp")



# -----------------------------
ptm = proc.time()[3]
aiksLiang = sum( apply(pAIKS_woWeight(Liang, score, c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
timeaiksLiang = proc.time()[3] - ptm
save(aiksLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', dd, '.RData'))
# save(aiksLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/b', nb, '/simAIKSLiang', dd, '.RData'))
# -----------------------------


# -----------------------------
# aiksLiang = sum( apply(pAIKS_woWeight(Liang, score, c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
# # save(aiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', dd, '.RData'))
# save(aiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/b', nb, '/simAIKSLiang', dd, '.RData'))
# 
# 
# cl = parallel::makeCluster(nprocs, type = "PSOCK")
# doParallel::registerDoParallel(cl)
# 
# nbs = unique(c(0, seq(nprocs, nb, by = nprocs), nb)); aiksBootLiang = c(); timeaiksLiang = 0
# # load(paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', dd, '.RData'))
# # nbs = unique(c(seq(length(aiksBootLiang), nb, by = nprocs), nb))
# 
# for(i in 2:length(nbs)){
#   ptm = proc.time()[3]
#   dummy = foreach(m = 1:(nbs[i]-nbs[i-1]), .combine = 'c', .packages = 'Rcpp', .noexport = c('AIKS', 'pAIKS')) %dopar% {
#     sourceCpp("RcppFtns.cpp")
#     qn = rmultinom(1, niter, rep(1/niter, niter))/niter - 1/niter
#     sum( apply(pAIKS(Liang, score, qn, c, beta, niter, 1), 2, sum) )
#   }
#   timeaiksLiang = timeaiksLiang + proc.time()[3] - ptm
#   aiksBootLiang = c(aiksBootLiang, dummy)
#   save(aiksLiang, aiksBootLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', dd, '.RData'))
#   # save(aiksLiang, aiksBootLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/b', nb, '/simAIKSLiang', dd, '.RData'))
# }
# -----------------------------
