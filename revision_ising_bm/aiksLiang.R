rm(list=ls())
library(Rcpp); library(RcppArmadillo)
library(doParallel); library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#================================================================
# Call data and functions
#================================================================
Nin = 1
# Nin = 2
# Nin = 3
# Nin = 4
# Nin = 5
# Nin = 6
# Nin = 7

load(paste('ACDAIKS/Liang/simAppxLiang', Nin, '.RData', sep = ''))




#================================================================
# AIKS
#================================================================
nchains = length(appx)
score = foreach(i = 1:nchains, .combine = c) %do% { appx[[i]][,1] }
nrep = foreach(i = 1:nchains, .combine = c) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j] == Liang[,i])) }
score = rep(score, nrep)
Liang = foreach(i = 1:nchains, .combine = c) %do% { Liang[,i] }

nn = 120000
# nn = 100000
# nn = 50000
indice = seq(20, nn, by = 20) # thinning
score = score[indice]
Liang = Liang[indice]
niter = length(Liang)

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

# -----------------------------
# nprocs = 19
# 
# ptm = proc.time()[3]
# aiksLiang = sum( apply(pIsingAIKS_woWeight(matrix(Liang, ncol = 1), matrix(score, ncol = 1), c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )
# timeaiksLiang = proc.time()[3] - ptm
# 
# save(aiksLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', Nin, '.RData'))
# -----------------------------


# -----------------------------
nprocs = 19
# nprocs = 7
mp_type = "PSOCK"
cl = parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

aiksLiang = sum( apply(pIsingAIKS_woWeight(matrix(Liang, ncol = 1), matrix(score, ncol = 1), c, beta, niter, nprocs+1), 2, sum)/(niter*(niter-1)) )

ptm = proc.time()[3]
aiksBootLiang = foreach(m = 1:nb, .combine = 'c', .packages = 'Rcpp', .noexport = c('IsingAIKS', 'pIsingAIKS')) %dopar% {
  sourceCpp("RcppFtns.cpp")
  qn = rmultinom(1, niter, rep(1/niter, niter))/niter - 1/niter
  sum( apply(pIsingAIKS(matrix(Liang, ncol = 1), matrix(score, ncol = 1), qn, c, beta, niter, 1), 2, sum) )
}
timeaiksLiang = proc.time()[3] - ptm

save(aiksLiang, aiksBootLiang, timeaiksLiang, file = paste0('ACDAIKS/Liang/n', nn, '/simAIKSLiang', Nin, '.RData'))
# -----------------------------






