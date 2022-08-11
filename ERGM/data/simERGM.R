rm(list=ls())
library(coda)
library(ergm)
library(Rcpp)
library(RcppArmadillo)
#library(MASS)
library(snow)
library(doParallel)
library(foreach)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


#======================================================================
# Call functions 
#======================================================================
sourceCpp("RcppFtns.cpp")



set.seed(1)
#======================================================================
# SIMULATION DTAT
#======================================================================
n = 30
data = network(n, density = 0.3, directed = FALSE) 
X = data[,]
formula = data ~ edges + gwesp(0.25,fixed=TRUE)
stat = summary(formula)
stat; Summary(X)
m = ergm(formula, estimate="MPLE")
hat = m$coef
summary(m)
COV = solve(-m$hessian)


### set true parameter values
Liang = ergmDMH(X, COV, matrix(hat, ncol = 2), 15000, 20)[-1,]
trth = round(colMeans(Liang[-(1:5000),]), 2)
X = Gibbs2(X, trth, 10000)
stat = Summary(X)
data = network(X, directed = FALSE)
# data = simulate(~ edges + gwesp(0.25,fixed=TRUE), nsim = 1, coef = trth, basis = net, 
#                 control = control.simulate.formula(MCMC.burnin = 10000))
# X = data[,]



### Parameter estimation via MLE or MPLE
formula = data ~ edges + gwesp(0.25,fixed=TRUE)
stat = summary(formula)
stat; Summary(X)
m = ergm(formula, estimate="MPLE")
hat = m$coef
summary(m)
COV = solve(-m$hessian)

save(X, data, stat, hat, COV, trth, file = "data/simERGM.RData")


