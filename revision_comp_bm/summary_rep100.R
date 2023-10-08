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
source("RFtns.R")


al = 0.01
# n = 400000
n = 300000
N = 200000
ncores = 20

summax = 8
# summax = 9
# summax = 10


# =========================================================
# ACD
# =========================================================
filenames = list.files(path = paste0("rep100/Trunc", summax), pattern = ".RData", all.files = TRUE, full.names = TRUE)

# ACD_rep100 = rep(0, 100); indices = 1:100
load(paste0('rep100/Trunc', summax, '/bids2ACDTrunc', summax, '.RData'))
indices = which(ACD_rep100 == 0)

# excl = c(58, 46, 59, 49, 50, 48, 51, 55, 52, 56, 57, 47, 54, 53, 60)
# excl = c(49, 56)
# indices = indices[!(indices %in% excl)]

for(i in indices){
  fname = paste0("rep100/Trunc", summax, "/rep", i, "_bids2AppxTrunc", summax, ".RData")
  if(fname %in% filenames){
    load(fname)

    p = ncol(Trunc)

    nrep = sapply(1:nth, function(j) sum(th[j,1]==Trunc[,1]))
    dhat = appx[, (p + 1):(p + p*(p+1)/2)]
    dhat = sapply(1:(p*(p+1)/2), function(j) rep(dhat[,j], nrep))

    ACD_rep100[i] = ACD_bm(dhat[1:n,], N)

    save(ACD_rep100, file = paste0('rep100/Trunc', summax, '/bids2ACDTrunc', summax, '.RData'))
  }
}



# =========================================================
# AIKS
# =========================================================
c = 1
beta = -1/2
nb = 1000

filenames = list.files(path = paste0("rep100/Trunc", summax), pattern = ".RData", all.files = TRUE, full.names = TRUE)

# AIKS_rep100 = rep(0, 100); indices = 1:100
load(paste0('rep100/Trunc', summax, '/bids2AIKSTrunc', summax, '.RData'))
indices = which(AIKS_rep100 == 0)

# excl = c(58, 46, 59, 49, 50, 48, 51, 55, 52, 56, 57, 47, 54, 53, 60)
# excl = c(49, 56)
# indices = indices[!(indices %in% excl)]

for(i in indices){
  fname = paste0("rep100/Trunc", summax, "/rep", i, "_bids2AppxTrunc", summax, ".RData")
  if(fname %in% filenames){
    load(fname)

    p = ncol(Trunc)
    score = appx[,1:p]
    nrep = sapply(1:nth, function(j) sum(th[j,1]==Trunc[,1]))
    score = sapply(1:p, function(j) rep(score[,j], nrep))
    
    indiceAIKS = seq(100, n, by = 100) # thinning
    score = score[indiceAIKS,]
    Trunc = Trunc[indiceAIKS,]
    niter = nrow(Trunc)

    sourceCpp("RcppFtns.cpp")
    AIKS_rep100[i] = sum( apply(pAIKS_woWeight(Trunc, score, c, beta, niter, ncores), 2, sum)/(niter*(niter-1)) )

    save(AIKS_rep100, file = paste0('rep100/Trunc', summax, '/bids2AIKSTrunc', summax, '.RData'))
  }
}




