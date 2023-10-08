rm(list=ls())
library(fields); library(lattice)
library(Rcpp); library(RcppArmadillo)
library(coda); library(xtable)
library(doParallel); library(foreach)
library(tidyverse); library(egg)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

source("RFtns.R")

#========================================================================
# Choose N
#========================================================================
# load("data/sim020Ising.RData")
# 
# N = 500000
# burn = 1000
# 
# sourceCpp("RcppFtns.cpp")
# Sx = Energy(X)
# theta = MPLE
# Sy = GibbStat(X, theta, N+burn)[-(1:burn)]
# 
# knots = seq(500, N, by = 100)
# appxN = foreach(j = knots, .combine = rbind) %do% {
#   uhat = f_uhat_b(Sx, Sy[1:j])
#   c(j, uhat, f_Hhat_b(Sx, Sy[1:j]) + uhat * uhat)
# }
# 
# p = 1
# datN = data.frame()
# for(i in 1:p){
#   datN = rbind(datN, data.frame(Name = paste0('u', i), Knot = appxN[,1], Value = appxN[,i+1]))
# }
# for(i in 1:(p*(p+1)/2)){
#   datN = rbind(datN, data.frame(Name = paste0('d', i), Knot = appxN[,1], Value = appxN[,i+p+1]))
# }


# save(datN, file = 'chooseNn/choose_N.RData')



#========================================================================
# Choose n
#========================================================================
Nin = 2
N = 200000
load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))

nchains = length(appx)
D = foreach(i = 1:nchains) %do% { appx[[i]][,2] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j]==Liang[, i])) }
D = foreach(i = 1:nchains, .combine = c) %do% { rep(D[[i]], nrep[[i]]) }

knots = seq(500, length(D), by = 100)
appxn = foreach(j = knots, .combine = rbind) %do% {
  c(j, f_Vhat_bm(D[1:j], N))
}

datn = data.frame()
for(i in 1:(ncol(appxn)-1)){
  datn = rbind(datn, data.frame(Name = paste0('V', i), Knot = appxn[,1], Value = appxn[,i+1]))
}

save(datn, file = 'chooseNn/choose_n.RData')



# ========================================================================
# plot
# ========================================================================
load('chooseNn/choose_N.RData')
load('chooseNn/choose_n.RData')


plot_aiksN = datN %>%
  filter(Name == 'u1', Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  labs(x = 'N', y = 'Approximation', title = expression('(a) '*u(theta)))
# plot_aiksN

plot_acdN = datN %>%
  filter(Name == 'd1', Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  labs(x = 'N', y = 'Approximation', title = expression('(b) '*d(theta)))
# plot_acdN

plot_acdn = datn %>%
  filter(Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  labs(x = 'n', y = 'Approximation', title = expression('(c) '*Sigma))
# plot_acdn


plot.final = ggarrange(plot_aiksN, plot_acdN, plot_acdn, nrow = 1)

ggsave(plot = plot.final, width = 9.2, height = 2.2, filename = 'figures_paper/isingChoosenN.eps')



