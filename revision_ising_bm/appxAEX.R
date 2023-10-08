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


#========================================================================
# Call a dataset, saved images, and functions
#========================================================================
dd = 100
# dd = 400
# dd = 450
# dd = 490
# dd = 500
# dd = 510

# appx = list()
# start = 1; timeappx = 0
load(paste0('ACDAIKS/AEX/simAppxAEX', dd, '.RData'))
start = length(appx) + 1


load("data/sim020Ising.RData")
load(paste0('ACDAIKS/AEX/simAEX', dd, '.RData'))



#========================================================================
# biased but consistent approximation for ACD
#========================================================================
niter = length(AEX)
# indice = burnin + 1:100000 # without thinning
indice = burnin + 1:200000 # without thinning
AEX = AEX[indice]

# nsets = 5
nsets = 10
AEX = matrix(AEX, ncol = nsets)
th = sapply(1:nsets, function(i) unique(AEX[,i]), simplify = F)
nth = sapply(1:nsets, function(i) length(th[[i]]))


### approximate
N = 200000
burn = 1000

sourceCpp("RcppFtns.cpp")
Sx = Energy(X)


nprocs = 19
# nprocs = 7
mp_type = "PSOCK"
cl = parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)


for(j in start:nsets){
  
  ptm = proc.time()[3]
  appx[[j]] = foreach(i = 1:(nth[j]), .combine = 'rbind', .packages = "Rcpp", .noexport = c("Energe", "GibbStat")) %dopar% {
    source("RFtns.R")
    sourceCpp("RcppFtns.cpp")
    
    Sy = GibbStat(X, th[[j]][i], N+burn)[-(1:burn)]
    c(f_uhat_b(Sx, Sy), f_dhat_b(Sx, Sy))
  }
  timeappx = timeappx + proc.time()[3] - ptm
  
  save(AEX, th, nth, Sx, appx, timeappx, file = paste0('ACDAIKS/AEX/simAppxAEX', dd, '.RData'))
}

