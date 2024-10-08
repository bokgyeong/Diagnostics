rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(snow)
library(doParallel)
library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")



# ==============================================================================
# Exchange algorithm
# ==============================================================================
load('data/bids2.RData')

betaInd = 1:p
beta_initial = rnorm(p)
nu = 1.75214 # MLE of nu by 'mpcmp'

b = rep(1, p)
block = sapply(1:length(unique(b)), function(i) which(b == i) - 1, simplify = F)
k = length(block)

sigma2 = rep(1^2, k)
COV = sapply(1:length(block), function(i) diag( length(block[[i]]) ), simplify = F)


### run real
thin = 1
burn = 10000
# niters = 100000 * thin + burn
niters = 500000 * thin + burn
updateCOV = TRUE
adaptInterval = 200
adaptFactorExponent = 0.8
adapIter = 1


# Murray = c(); rtime = 0; preniters = 0; Accprob = 0

load(paste0('ACDAIKS/Murray_knu/bids2Murray.RData'))
preniters = nrow(Murray)

sourceCpp('RcppFtns.cpp')

ptm = proc.time()[3]
dummy = exactCOMPreg_knu(outer = niters - preniters, y = y, X = X, beta = beta_initial, nu = nu,
                         block = block, sigma2 = sigma2, COV = COV, updateCOV = updateCOV, 
                         adaptInterval = adaptInterval, adaptFactorExponent = adaptFactorExponent,
                         adapIter = adapIter, thin = thin)
rtime = rtime + proc.time()[3] - ptm

Murray = rbind(Murray, dummy$Sample)
Accprob = ( Accprob * preniters + colSums(dummy$Accprob) ) / niters

nSamples = nrow(dummy$Sample)
beta_initial = dummy$Sample[nSamples, betaInd]
sigma2 = dummy$sigma2
adapIter = dummy$adapIter
COV = dummy$COV

save(Murray, Accprob, rtime, thin, block, burn,
     sigma2, adapIter, COV, beta_initial, 
     updateCOV, adaptInterval, adaptFactorExponent, nu,
     file = paste0('ACDAIKS/Murray_knu/bids2Murray.RData'))

