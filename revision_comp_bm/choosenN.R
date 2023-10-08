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
# load('data/bids2.RData')
# load(paste0('ACD/Trunc_knu/bids2Trunc', 10, '.RData'))
# 
# N = 500000
# 
# sourceCpp("RcppFtns.cpp")
# theta = Trunc[nrow(Trunc),]
# Sx = t(X) %*% y
# nu = 1.75214
# Sy = rCOMP2_parallel_knu(X, theta, nu, N, 1)
# 
# knots = seq(500, N, by = 100)
# appxN = foreach(j = knots, .combine = rbind) %do% {
#   Uhat = f_Uhat_knu(Sx, Sy[1:j,], nu)
#   dhat = f_dhat_knu(Sx, Sy[1:j,], nu, Uhat)
#   c(j, Uhat, dhat)
# }
# 
# p = 10
# datN = data.frame()
# for(i in 1:p){
#   datN = rbind(datN, data.frame(Name = paste0('u', i), Knot = appxN[,1], Value = appxN[,i+1]))
# }
# for(i in 1:(p*(p+1)/2)){
#   datN = rbind(datN, data.frame(Name = paste0('d', i), Knot = appxN[,1], Value = appxN[,i+p+1]))
# }
# 
# save(datN, file = 'choosenN/choose_N.RData')



#========================================================================
# Choose n
#========================================================================
Nin = 8
load(paste0('ACDAIKS/Trunc_knu/bids2AppxTrunc', Nin, '.RData'))

N = 200000
nsets = length(appx)
p = ncol(Trunc[[1]])

nrep = foreach(i = 1:nsets) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Trunc[[i]][,1] & th[[i]][j,p]==Trunc[[i]][,p])) }
dhat = foreach(i = 1:nsets) %do% { appx[[i]][, (p + 1):(p + p*(p+1)/2)] }
dhat = foreach(i = 1:nsets, .combine = rbind) %do% { sapply(1:(p*(p+1)/2), function(j) rep(dhat[[i]][,j], nrep[[i]])) }

knots = seq(500, nrow(dhat), by = 100)
appxn = foreach(j = knots, .combine = rbind) %do% {
  Vhat = f_Vhat_bm(dhat[1:j,], N)
  c(j, Vhat[upper.tri(Vhat, diag = T)])
}

datn = data.frame()
# for(i in 1:(ncol(appxn)-1)){
for(i in 1:10){
  datn = rbind(datn, data.frame(Name = paste0('V', i), Knot = appxn[,1], Value = appxn[,i+1]))
}

save(datn, file = 'choosenN/choose_n.RData')



# ========================================================================
# plot
# ========================================================================
load('choosenN/choose_N.RData')
load('choosenN/choose_n.RData')

datn$Name = factor(datn$Name, levels = paste0('V', 1:10), labels = c(paste0('S1', 1:10)))


plot_aiksN = datN %>%
  filter(Name %in% paste0('u', 1:3), Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  facet_wrap(~Name, nrow = 1, scales = 'free') +
  labs(x = 'N', y = 'Approximation', title = expression('(a) '*u(theta)))
# plot_aiksN

plot_acdN = datN %>%
  filter(Name %in% paste0('d', 1:3), Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  facet_wrap(~Name, nrow = 1, scales = 'free') +
  labs(x = 'N', y = 'Approximation', title = expression('(b) '*d(theta)))
# plot_acdN

plot_acdn = datn %>%
  filter(Name %in% paste0('S1', 1:3), Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  facet_wrap(~Name, nrow = 1, scales = 'free') +
  labs(x = 'n', y = 'Approximation', title = expression('(c) '*Sigma))
# plot_acdn


plot.final = ggarrange(plot_aiksN, plot_acdN, plot_acdn, ncol = 1)

ggsave(plot = plot.final, width = 9.2, height = 7, filename = 'figures_paper/compChoosenN.eps')



