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
dd = 10
# dd = 20
# dd = 50
# dd = 100
# dd = 200
# dd = 400

# appx = list()
# start = 1; timeappx = 0
load(paste0('ACDAIKS/Atchade/simAppxAtchade', dd, '.RData'))
start = length(appx) + 1

load("data/sim020Ising.RData")
load(paste0('ACDAIKS/Atchade/simAtchade', dd, '.RData'))



#========================================================================
# biased but consistent approximation for ACD
#========================================================================
niter = length(Atchade)
# indice = burnin + 1:100000 # without thinning
indice = burnin + 1:200000 # without thinning
Atchade = Atchade[indice]

# nchains = 5
nchains = 10
Atchade = matrix(Atchade, ncol = nchains)
th = sapply(1:nchains, function(i) unique(Atchade[,i]), simplify = F)
nth = sapply(1:nchains, function(i) length(th[[i]]))


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



for(j in start:nchains){
  
  ptm = proc.time()[3]
  appx[[j]] = foreach(i = 1:(nth[j]), .combine = 'rbind', .packages = "Rcpp", .noexport = c("Energe", "GibbStat")) %dopar% {
    source("RFtns.R")
    sourceCpp("RcppFtns.cpp")
    
    Sy = GibbStat(X, th[[j]][i], N+burn)[-(1:burn)]
    c(f_uhat_b(Sx, Sy), f_dhat_b(Sx, Sy))
  }
  timeappx = timeappx + proc.time()[3] - ptm
  
  save(Atchade, th, nth, Sx, appx, timeappx, file = paste0('ACDAIKS/Atchade/simAppxAtchade', dd, '.RData'))
}

