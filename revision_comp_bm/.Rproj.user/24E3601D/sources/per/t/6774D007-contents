rm(list=ls())
library(Rcpp); library(RcppArmadillo)
library(doParallel); library(foreach)
library(tidyverse); library(egg); library(gridExtra)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
source("RFtns.R")


#========================================================================
# load data
#========================================================================
load('data/bids2.RData')

summax = 10
load(paste0('ACDAIKS/Trunc_knu/bids2Trunc', summax, '.RData'))


rtime = list()
nrepcost = 100
npostACD = 73000
npostAIKS = 4000
N = 200000
ncores = 20


#========================================================================
# approximation
#========================================================================
sourceCpp("RcppFtns.cpp")


### simulation --------
ptm = proc.time()[3]
Sy = rCOMP2_parallel_knu(X, Trunc[nrow(Trunc),], nu, N, 1)
etm = proc.time()[3] - ptm
rtime$simulation = etm / N


data.cost = data.frame(n = npostACD, nprime = npostAIKS, sim = rtime$simulation)


#========================================================================
# ACD
#========================================================================
summax = 10
load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))

nsets = length(appx)
p = ncol(Trunc[[1]])

nrep = foreach(i = 1:nsets) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Trunc[[i]][,1] & th[[i]][j,p]==Trunc[[i]][,p])) }
dhat = foreach(i = 1:nsets) %do% { appx[[i]][, (p + 1):(p + p*(p+1)/2)] }
dhat = foreach(i = 1:nsets, .combine = rbind) %do% { sapply(1:(p*(p+1)/2), function(j) rep(dhat[[i]][,j], nrep[[i]])) }

source('RFtns.R')

dummy = c(); k = 1
for(n in unique(data.cost$n)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    acdTrunc = ACD_bm(dhat[1:n,], N)
    etm = etm + proc.time()[3] - ptm
  }
  dummy[k] = etm / nrepcost
  k = k + 1
}
data.cost$acd = dummy



#========================================================================
# AIKS
#========================================================================
summax = 10
load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', summax, '.RData'))


p = ncol(Trunc[[1]])
nchains = length(appx)
score = foreach(i = 1:nchains) %do% { appx[[i]][,1:p] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1] == Trunc[[i]][,1])) }
score = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:p, function(j) rep(score[[i]][,j], nrep[[i]])) }
Trunc = foreach(i = 1:nchains, .combine = rbind) %do% { Trunc[[i]] }

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

dummy = c(); k = 1
for(n in unique(data.cost$nprime)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    aiksTrunc = sum( apply(pAIKS_woWeight(Trunc[1:n,], score[1:n,], c, beta, n, ncores), 2, sum)/(n*(n-1)) )
    etm = etm + proc.time()[3] - ptm
  }
  dummy[k] = etm / nrepcost
  k = k + 1
}
save(data.cost, dummy, file = 'cost/cost_diag.RData')

data.cost$aiks = mean(dummy)
save(data.cost, file = 'cost/cost_diag.RData')






# ==============================================================================
# plot
# ==============================================================================
load('cost/cost_diag.RData')

nb = 1000
ncores = 20
N = 200000
data.cost.total = data.frame()
for(j in N){
  dummy = data.cost %>%
    add_column(N = j) %>%
    mutate(simACD = paste0(round(sim*n*N/ncores/60/60, 2), ' hr'),
           simAIKS = paste0(round(sim*nprime*N/ncores/60, 2), ' min'),
           ACD = paste0(round(acd, 4), ' sec'), AIKS = paste0(round(aiks/60, 2), ' min')) %>%
    select(n, N, simACD, simAIKS, ACD, AIKS)
  data.cost.total = rbind(data.cost.total, dummy)
}
data.cost.total

