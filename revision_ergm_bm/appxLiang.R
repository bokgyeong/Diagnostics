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



#========================================================================
# call data and functions
#========================================================================
# Nin = 1
# Nin = 2
# Nin = 3
# Nin = 4
# Nin = 5
# Nin = 6
Nin = 20


appx = list()
start = 1; timeappx = 0
# load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))
# start = length(appx) + 1


load(paste0('ACDAIKS/Liang/simLiang', Nin, '.RData'))




#========================================================================
# Decide lag
#========================================================================
# niter = nrow(Liang)
# 
# res = acf(Liang[(burnin+1):niter,])$acf
# ACF = cbind(res[,,1], res[,,2])
# 
# ACF = ACF[-c(1),]
# Lags = acf(Liang[(burnin+1):niter,])$lag
# Lags = cbind(Lags[-c(1),,1],Lags[-c(1),,2])
# Names = c('theta 1', 'theta 2 & theta 1', 'theta 1 & theta 2', 'theta 2')
# Threshs = c(0.1, -0.1, -0.1, 0.1)
# 
# dat.thin = data.frame()
# for(i in 1:ncol(ACF)){
#   dat = data.frame(Name = Names[i], Lag = Lags[,i], ACF = ACF[,i], Thresh = Threshs[i])
#   dat.thin = rbind(dat.thin, dat)
# }
# 
# library(tidyverse)
# plot.thin = dat.thin %>% ggplot(aes(Lag, ACF, group = 1)) +
#   geom_point() +
#   geom_line() +
#   geom_hline(aes(yintercept = Thresh), linetype = 'dashed', color = 'blue') +
#   facet_wrap(~ Name, scales = 'free') +
#   labs(x = 'Lag', y = 'Average ACF (10 chains)')
# plot.thin
# 
# ggsave(filename = 'figures/simLiangACF.eps', plot = plot.thin, device = cairo_ps, width = 8, height = 4.5)

# 20 looks good


#========================================================================
# biased but consistent approximation
#========================================================================
niter = nrow(Liang)
# indice = burnin + 1:100000 # without thinning
# indice = burnin + 1:500000 # without thinning
indice = burnin + 1:250000 # without thinning
Liang =  Liang[indice,]
niter = nrow(Liang)

# nsets = 10
nsets = 25
# nsets = 50
Liang = sapply(1:nsets, function(i) Liang[((i-1) * (niter/nsets) + 1):(i * niter/nsets),], simplify = F)
th = sapply(1:nsets, function(i) unique(Liang[[i]]), simplify = F)
nth = sapply(1:nsets, function(i) nrow(th[[i]]))


### simulation by Gibbs sampler
N = 200000
burn = 1000

nprocs = 19
# nprocs = 7
mp_type = "PSOCK"
cl = parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)


for(j in start:nsets){
  
  ptm = proc.time()[3]
  appx[[j]] = foreach(i = 1:(nth[j]), .combine = 'rbind', .packages = "Rcpp") %dopar% {
    source("RFtns.R")
    sourceCpp("RcppFtns.cpp")
    Sy = Gibbs3(X, th[[j]][i,], N+burn)[-(1:burn),]
    c(f_Uhat(stat, Sy), f_dhat(stat, Sy))
  }
  timeappx = timeappx + proc.time()[3] - ptm
  
  save(Liang, th, nth, stat, appx, timeappx, file = paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))
}



# ptm = proc.time()[3]
# Sy = Gibbs3(X, th[[1]][1,], 1000)
# rtime = proc.time()[3] - ptm
# rtime/1000


### Choose N -------------------------------------------------------------------

# N = 500000
# burn = 1000
# 
# sourceCpp("RcppFtns.cpp")
# Sy = Gibbs3(X, Liang[1,], N+burn)[-(1:burn),]
# 
# source("RFtns.R")
# knots = seq(100, N, by = 100)
# appxN = foreach(j = knots, .combine = rbind) %do% {
#   c(j, f_Hhat(stat, Sy[1:j,]) + f_Jhat(stat, Sy[1:j,]))
# }
# 
# datN = data.frame()
# for(i in 1:(ncol(appxN)-1)){
#   datN = rbind(datN, data.frame(Name = paste0('d', i), Knot = appxN[,1], Value = appxN[,i+1]))
# }
# 
# library(tidyverse)
# plot_N = datN %>%
#   filter(Knot > 10000) %>%
#   ggplot(aes(Knot, Value)) +
#   geom_line() +
#   facet_wrap(~ Name, scales = 'free') +
#   labs(x = 'N', y = 'Approximation')
# plot_N
# 
# # save(datN, file = 'figures/choose_N.RData')

# N = 200,000 looks good
# N = 300,000 looks better
### ----------------------------------------------------------------------------


