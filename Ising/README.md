Ising Model
========


Install required R packages and call R and Rcpp functions
---------------------------

``` r
# Install required R packages
library(fields)
library(lattice)
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(xtable)
library(doParallel)
library(foreach)
library(batchmeans)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

# Call R and Rcpp functions
source("RFtns.R")
sourceCpp("RcppFtns.cpp")
```

Set the R working directory to the location of this README file.


Simulate data from the Ising Model
---------------------------

``` r
# Simulate N by N lattice from the Ising model
set.seed(1)
N = 30
X = ProppWilson(N, N, 0.20) # Perfect sampling from the Ising model with dependence parameter of 0.20

# Maximum pseudolikelihood estimate for initial value
z = optim(0.1, pseudoLik, control = list(fnscale = -1), method = "BFGS")
MPLE = z$par
```


Generate posterior sample via double Metropolis-Hastings (DMH) algorithm
---------------------------

``` r
Nin = 1 # Length of inner sampler of DMH
burnin = 1000
Nout = 120,000 + burnin # Length of outer sampler of DMH

ptm = proc.time()[3]
Liang = IsingDMH(Nout, Nin, MPLE, 0.1, X) # Standard deviation of normal proposal = 0.1
timeLiang = proc.time()[3] - ptm
```


Approximate intractable terms in our diagnostics
---------------------------

``` r
th = unique(Liang) # unique posterior sample points
nth = length(th)

N = 200000 # Number of auxilairy samples for Monte Carlo approximation
burn = 1000
Sx = Energy(X) # Sufficient statistic

# Set up parallelization for 20 cores
cl = parallel::makeCluster(19, type = "PSOCK")
doParallel::registerDoParallel(cl)

# Approximate intractable terms in parellel
ptm = proc.time()[3]
appx = foreach(i = 1:nth, .combine = 'rbind', .packages = "Rcpp", .noexport = c("Energe", "GibbStat")) %dopar% {
  source("RFtns.R")
  sourceCpp("RcppFtns.cpp")
  
  Sy = GibbStat(X, th, N+burn)[-(1:burn)]
  c(f_uhat_b(Sx, Sy), f_dhat_b(Sx, Sy))
}
timeappx = proc.time()[3] - ptm
```


Compute approximate curvature diagnostic (ACD)
---------------------------

``` r
# Approximate difference between hassian and outerproject of score at unique posterior sample points
D = appx[,2] 

# Approximate difference between hassian and outerproject of score at entire posterior sample points
nrep = sapply(1:nth, function(j) sum(th[j] == Liang))
D = rep(D, nrep) 

# Compute ACD
ptm = proc.time()[3]
acdLiang = ACD_bm(D)
timeacdLiang = proc.time()[3] - ptm

# Threshold of ACD
qchisq(0.99, 1)
```


Compute approximate inverse multiquadric kernel Stein discrepancy (AIKS)
---------------------------

``` r
# Approximate score function at unique posterior sample points
score = appx[,1]

# Approximate score function at entire posterior sample points
nrep = sapply(1:nth, function(j) sum(th[j] == Liang))
score = rep(score, nrep)

# Thinning
indice = seq(20, length(Liang), by = 20) 
score = score[indice]
Liang = Liang[indice]
niter = length(Liang)

# Tuning parameters of AIKS
c = 1
beta = -1/2

# Compute AIKS in parallel for 20 cores
ptm = proc.time()[3]
aiksLiang = sum( apply(pIsingAIKS_woWeight(matrix(Liang, ncol = 1), matrix(score, ncol = 1), c, beta, niter, 20), 2, sum) / (niter * (niter-1)) )
timeaiksLiang = proc.time()[3] - ptm

# Compute AIKS threshold in parallel for 20 cores
nb = 1000 # Bootstrap sample size

cl = parallel::makeCluster(19, type = "PSOCK")
doParallel::registerDoParallel(cl)

ptm = proc.time()[3]
aiksBootLiang = foreach(m = 1:nb, .combine = 'c', .packages = 'Rcpp', .noexport = c('IsingAIKS', 'pIsingAIKS')) %dopar% {
  sourceCpp("RcppFtns.cpp")
  qn = rmultinom(1, niter, rep(1/niter, niter))/niter - 1/niter
  sum( apply(pIsingAIKS(matrix(Liang, ncol = 1), matrix(score, ncol = 1), qn, c, beta, niter, 1), 2, sum) )
}
timeaiksBootLiang = proc.time()[3] - ptm

quantile(aiksBootLiang, 0.99) # Threshold of AIKS

```

