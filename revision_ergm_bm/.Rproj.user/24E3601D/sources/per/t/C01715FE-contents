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
# load("data/simERGM.RData")
# 
# N = 500000
# burn = 1000
# 
# sourceCpp("RcppFtns.cpp")
# theta = hat
# Sy = Gibbs3(X, theta, N+burn)[-(1:burn),]
# 
# knots = seq(500, N, by = 100)
# appxN = foreach(j = knots, .combine = rbind) %do% {
#   c(j, f_Uhat(stat, Sy[1:j,]), f_Hhat(stat, Sy[1:j,]) + f_Jhat(stat, Sy[1:j,]))
# }
# 
# p = 2
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
Nin = 2
load(paste0('ACDAIKS/Liang/simAppxLiang', Nin, '.RData'))

N = 200000
nchains = length(appx)
D = foreach(i = 1:nchains) %do% { appx[[i]][,3:5] }
nrep = foreach(i = 1:nchains) %do% { sapply(1:nth[i], function(j) sum(th[[i]][j,1]==Liang[[i]][,1])) }
D = foreach(i = 1:nchains, .combine = rbind) %do% { sapply(1:3, function(j) rep(D[[i]][,j], nrep[[i]])) }

knots = seq(500, nrow(D), by = 100)
appxn = foreach(j = knots, .combine = rbind) %do% {
  Vhat = f_Vhat_bm(D[1:j,], N)
  c(j, Vhat[upper.tri(Vhat, diag = T)])
}

datn = data.frame()
for(i in 1:(ncol(appxn)-1)){
  datn = rbind(datn, data.frame(Name = paste0('V', i), Knot = appxn[,1], Value = appxn[,i+1]))
}

save(datn, file = 'choosenN/choose_n.RData')



# ========================================================================
# plot
# ========================================================================
load('choosenN/choose_N.RData')
load('choosenN/choose_n.RData')

datn$Name = factor(datn$Name, levels = paste0('V', 1:3), labels = c('S11', 'S12', 'S22'))

plot_aiksN = datN %>%
  filter(Name %in% paste0('u', 1:2), Knot > 1000) %>%
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
  filter(Name %in% c('S11', 'S12', 'S22'), Knot > 1000) %>%
  ggplot(aes(Knot, Value)) +
  geom_line() +
  facet_wrap(~Name, nrow = 1, scales = 'free') +
  labs(x = 'n', y = 'Approximation', title = expression('(c) '*Sigma))
# plot_acdn


plot.blank = ggplot() + theme_void()
plot_aiksN_blank = ggarrange(plot_aiksN, plot.blank, nrow = 1, widths = c(1.92, 1))

plot.final = grid.arrange(plot_aiksN_blank, plot_acdN, plot_acdn, ncol = 1)

ggsave(plot = plot.final, width = 9.2, height = 7, filename = 'figures_paper/ergmChoosenN.eps')



