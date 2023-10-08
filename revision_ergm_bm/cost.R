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
load("data/simERGM.RData")

Nin = 1
load(paste0('ACDAIKS/Liang/simLiang', Nin, '.RData'))


rtime = list()
nrepcost = 100
npostACD = 106000
npostAIKS = 10000
N = 200000
ncores = 20

#========================================================================
# approximation
#========================================================================
sourceCpp("RcppFtns.cpp")


### simulation --------
ptm = proc.time()[3]
Sy = Gibbs3(X, Liang[nrow(Liang),], N)
etm = proc.time()[3] - ptm
rtime$simulation = etm / N


data.cost = data.frame(n = npostACD, nprime = npostAIKS, sim = rtime$simulation)


#========================================================================
# ACD
#========================================================================
Nin = 1
load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))
nchains = length(appx)

D = foreach(i = 1:nchains) %do% { appx[[i]][,3:5] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Liang[[i]][,1])) }
D = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:3, function(j) rep(D[[i]][,j], nrep[[i]])) }

source('RFtns.R')

dummy = c(); k = 1
for(n in unique(data.cost$n)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    acdLiang = ACD_bm(D[1:n,])
    etm = etm + proc.time()[3] - ptm
  }
  dummy[k] = etm / nrepcost
  k = k + 1
}
data.cost$acd = dummy



#========================================================================
# AIKS
#========================================================================
Nin = 1
load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))

p = ncol(Liang[[1]])
nchains = length(appx)
score = foreach(i = 1:nchains) %do% { appx[[i]][,1:p] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1] == Liang[[i]][,1])) }
score = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:p, function(j) rep(score[[i]][,j], nrep[[i]])) }
Liang = foreach(i = 1:nchains, .combine = rbind) %do% { Liang[[i]] }

nn = 250000
indice = seq(25, nn, by = 25) # thinning
score = score[indice,]
Liang = Liang[indice,]
niter = nrow(Liang)

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

dummy = c(); k = 1
for(n in unique(data.cost$nprime)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    aiksLiang = sum( apply(pAIKS_woWeight(Liang[1:n,], score[1:n,], c, beta, n, ncores), 2, sum)/(n*(n-1)) )
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

# nb = 1000
# ncores = 20
# N = 200000
# data.cost.total = data.frame()
# for(j in N){
#   dummy = data.cost %>%
#     add_column(N = j) %>%
#     mutate(Simulation = sim * n * N / ncores / 60 / 60) %>%
#     mutate(ACD = Simulation + acd/60/60, AIKS = Simulation + aiks/60/60, AIKSthr = aiks*nb/60/60) %>%
#     select(n, N, Simulation, ACD, AIKS, AIKSthr)
#   data.cost.total = rbind(data.cost.total, dummy)
# }
# data.cost.total

nb = 1000
ncores = 20
N = 200000
data.cost.total = data.frame()
for(j in N){
  dummy = data.cost %>%
    add_column(N = j) %>%
    mutate(simACD = paste0(round(sim*n*N/ncores/60/60, 2), ' hr'),
           simAIKS = paste0(round(sim*nprime*N/ncores/60/60, 2), ' hr'),
           AIKSthred = paste0(round(aiks*nb/60/60, 2), ' hr'),
           ACD = paste0(round(acd, 4), ' sec'), AIKS = paste0(round(aiks/60, 2), ' min')) %>%
    select(n, N, Simulation, AIKSthred, ACD, AIKS)
  data.cost.total = rbind(data.cost.total, dummy)
}
data.cost.total


# data.cost.total = data.cost.total %>%  
#   select(n, N, ACD, AIKS) %>% 
#   pivot_longer(
#     cols = c('ACD', 'AIKS'), 
#     names_to = "Method", 
#     values_to = "Cost"
#   )
# 
# # data.cost.total$N = paste0('N = ', format(data.cost.total$N, scientific = F))
# # 
# # plot.cost = data.cost.total %>%
# #   ggplot(aes(x = n, y = Cost)) +
# #   geom_line(aes(color = Method, linetype = Method)) +
# #   geom_point(aes(color = Method, shape = Method), size = 2) +
# #   facet_wrap(~ N, nrow = 1) +
# #   labs(x = 'Posterior sample size, n', y = 'Computation cost (hour)')
# # 
# # ggsave(plot.cost, 
# #        filename = 'figures_slide/cost1.eps', width = 9, height = 2.5)
# 
# 
# data.cost.total$N = format(data.cost.total$N, scientific = F)
# 
# plot.cost2 = data.cost.total %>%
#   ggplot(aes(x = n, y = Cost)) +
#   geom_line(aes(color = N, linetype = N)) +
#   geom_point(aes(color = N, shape = N), size = 2) +
#   facet_wrap(~ Method, nrow = 1) +
#   labs(x = 'Posterior sample size, n', y = 'Computation cost (hour)')
# 
# ggsave(plot.cost2,
#        filename = 'figures_slide/cost2.eps', width = 5.3, height = 2.5)


