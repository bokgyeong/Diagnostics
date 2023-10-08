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
load("data/sim020Ising.RData")

Nin = 1
load(paste0('ACDAIKS/Liang/simLiang', Nin, '.RData'))


rtime = list()
nrepcost = 100
npostACD = 24000
npostAIKS = 6000
N = 200000
ncores = 20

#========================================================================
# approximation
#========================================================================
sourceCpp("RcppFtns.cpp")
Sx = Energy(X)


### simulation --------
ptm = proc.time()[3]
Sy = GibbStat(X, Liang[nrow(Liang),1], N)
etm = proc.time()[3] - ptm
rtime$simulation = etm / N


data.cost = data.frame(n = npostACD, nprime = npostAIKS, sim = rtime$simulation)


#========================================================================
# ACD
#========================================================================
Nin = 1
load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))
nchains = length(appx)

D = foreach(i = 1:nchains) %do% { appx[[i]][,2] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j]==Liang[, i])) }
D = foreach(i = 1:nchains, .combine = c) %do% { rep(D[[i]], nrep[[i]]) }

source('RFtns.R')

dummy = c(); k = 1
for(n in unique(data.cost$n)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    acdLiang = ACD_bm(D[1:n])
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

nchains = length(appx)
score = foreach(i = 1:nchains, .combine = c) %do% { appx[[i]][,1] }
nrep = foreach(i = 1:nchains, .combine = c) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j] == Liang[,i])) }
score = rep(score, nrep)
Liang = foreach(i = 1:nchains, .combine = c) %do% { Liang[,i] }

nn = 120000
indice = seq(20, nn, by = 20) # thinning
score = score[indice]
Liang = Liang[indice]
niter = length(Liang)

c = 1
beta = -1/2
nb = 1000
sourceCpp("RcppFtns.cpp")

dummy = c(); k = 1
for(n in unique(data.cost$nprime)){
  etm = 0
  for(i in 1:nrepcost){
    ptm = proc.time()[3]
    aiksLiang = sum( apply(pIsingAIKS_woWeight(matrix(Liang[1:n], ncol = 1), matrix(score[1:n], ncol = 1), c, beta, n, ncores), 2, sum)/(n*(n-1)) )
    etm = etm + proc.time()[3] - ptm
  }
  dummy[k] = etm / nrepcost
  k = k + 1
}

# cl = parallel::makeCluster(19, type = "PSOCK")
# doParallel::registerDoParallel(cl)
#
# dummy = foreach(n = unique(data.cost$n), .combine = 'cbind', .packages = 'Rcpp', .noexport = c('IsingAIKS_woWeight', 'pIsingAIKS_woWeight')) %:%
#   foreach(i = 1:nrepcost, .combine = 'c') %dopar% {
#     sourceCpp("RcppFtns.cpp")
#     ptm = proc.time()[3]
#     aiksLiang = sum( apply(pIsingAIKS_woWeight(matrix(Liang[1:n], ncol = 1), matrix(score[1:n], ncol = 1), c, beta, n, 1), 2, sum)/(n*(n-1)) )
#     proc.time()[3] - ptm
#   }

save(data.cost, dummy, file = 'cost/cost_diag.RData')

data.cost$aiks = mean(dummy)
save(data.cost, file = 'cost/cost_diag.RData')






# ==============================================================================
# plot
# ==============================================================================
load('cost/cost_diag.RData')

# nb = 1000
# ncores = 20
# N = c(50000, 100000, 200000, 300000)
# data.cost.total = data.frame()
# for(j in N){
#   dummy = data.cost %>%
#     add_column(N = j) %>%
#     mutate(Simulation = sim * n * N / ncores / 60 / 60) %>%
#     mutate(ACD = Simulation + acd/60/60, AIKS = Simulation + aiks*nb/ncores/60/60) %>%
#     select(n, N, Simulation, ACD, AIKS)
#   data.cost.total = rbind(data.cost.total, dummy)
# }

nb = 1000
ncores = 20
N = 200000
data.cost.total = data.frame()
for(j in N){
  dummy = data.cost %>%
    add_column(N = j) %>%
    mutate(simACD = paste0(round(sim*n*N/ncores/60/60, 2), ' hr'),
           simAIKS = paste0(round(sim*nprime*N/ncores/60/60, 2), ' hr'),
           AIKSthred = paste0(round(aiks*nb/ncores/60, 2), ' min'),
           ACD = paste0(round(acd, 4), ' sec'), AIKS = paste0(round(aiks, 2), ' sec')) %>%
    select(n, N, simACD, simAIKS, AIKSthred, ACD, AIKS)
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


