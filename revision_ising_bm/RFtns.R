### Bhattacharyya distance function -----------------------------
Bdist = function(a,b,n){
  p = density(a,n=2^15)
  q = density(b,n=2^15)
  sample = c()
  
  for(i in 1:n){
    u = runif(1)
    ind.p = which.min( abs( p$x - u ) )
    ind.q = which.min( abs( q$x - u ) ) 
    sample[i] = sqrt( p$y[ind.p]*q$y[ind.q] )}
  
  result = mean(sample)
  if( result < 1 ){ result = result }else{ 
    result = result - 2*(result - 1) }
  
  return(  -log(result)  )              
}

### sign corrected version of batch mean --------------------------------
bmsign = function(vals,sign,bs="sqroot",warn=FALSE){
  N <- length(vals)
  if (N<1000)
  {
    if (warn) # if warning
      cat("WARNING: too few samples (less than 1000)\n")
    if (N<10)
      return(NA)
  }
  
  if (bs=="sqroot") 
  {
    b <- floor(sqrt(N)) # batch size
    a <- floor(N/b) # number of batches
  }
  else
    if (bs=="cuberoot") 
    {
      b <- floor(N^(1/3)) # batch size
      a <- floor(N/b) # number of batches
    }
  else # batch size provided
  {
    stopifnot(is.numeric(bs))  
    b <- floor(bs) # batch size
    if (b > 1) # batch size valid
      a <- floor(N/b) # number of batches
    else
      stop("batch size invalid (bs=",bs,")")
  }
  
  Ys <- sapply(1:a,function(k) 
    return(  mean( vals[((k-1)*b+1):(k*b)]*sign[((k-1)*b+1):(k*b)] )/mean(sign[((k-1)*b+1):(k*b)])  ))
  
  muhat <- mean(vals*sign)/mean(sign)
  sigmahatsq <- b*sum((Ys-muhat)^2)/(a-1)
  
  bmse <- sqrt(sigmahatsq/N)
  
  return(list(est=muhat,se=bmse))
}


### log likelihood for MPLE -------------------------------------------
pseudoLik <- function(parameter){
  beta <- parameter
  work <- cbind(0, rbind(0, X, 0), 0); c <- 0
  for(i in 2:(N+1)){
    for(j in 2:(N+1)){
      p <- exp(2*beta*(work[i,j-1]+work[i,j+1]+work[i-1,j]+work[i+1,j]))
      p <- p/(1+p)
      c <- c + dbinom((work[i,j]+1)/2,1,p,log=TRUE)}}
  return(c)
}





### Particle generation for function emulation algorithm
partFDMH = function(X, stat, MPLE, nPar, num){ # d: no of particles
  parts = rep(0, nPar)
  Nout <- 20000
  Nin = 10
  
  FLiang = IsingFDMH(Nout, Nin, MPLE, 0.1, X)  # multiply 0.5 in c code
  # FLiang = IsingDMH(Nout, Nin, MPLE, 0.1, X)
  
  
  FLiang = FLiang[501:Nout]                               # burn in 500
  stand =  (FLiang-min(FLiang))/(max(FLiang)-min(FLiang)) # standardized
  stand = unique(stand)                                   # only take unique components
  dmat = rdist(stand)                                     # distance mat
  
  # choose auxiliary parameters through min max procedure
  ind = 1; A = 1; Ac = 2:length(stand)
  parts[1] = stand[ind]
  
  ind = which.max( dmat[,A] )
  A = c(A,ind)
  Ac = Ac[-which(Ac==ind)]
  parts[2] = stand[ind]
  
  
  for(i in 3:nPar){
    dummy = max( apply( dmat[,A] , 1, min )[Ac] )
    ind = which( dmat[,A] == dummy  )
    if(ind < dim(dmat)[1]){ ind = ind }else{ ind = ind-floor( ind/dim(dmat)[1] )*dim(dmat)[1] }
    A = c(A,ind)
    Ac = Ac[-which(Ac==ind)]
    parts[i] = stand[ind]
  }
  
  dist.parts = rdist(parts)  # distance matrix for parts (for standardized version)
  parts = (max(FLiang)-min(FLiang))*parts + min(FLiang)
  
  return(parts)
}


### Gaussian Process MCMC with ABC particle generation
# IS approximation with independent samples
IsingGPmcmc = function(Niter, d, N, X, stat, cycle = 100, LikEm = FALSE, num){ 
  p = 1
  Domain = c(0, 1)
  num.point = 3000           
  th = matrix(0,num.point,p)
  
  # Generate D Latin hypercube design points over the domain D_1
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2]-Domain[1])*th[,1]+Domain[1]
  
  # Simulate the auxiliary variable for each design point
  summat = pAuxStatGibb2(X, 100, th, num) # 100 cycle
  
  # Choose some of them based on the distance b/w simulated and observed summary statistics
  dist = sqrt( ( rep(stat, num.point) - summat[1:num.point,] )^2 )
  eps = quantile(dist, probs = 0.03)
  m =  mean(th[which(dist < eps),])
  S = var( th[which(dist < eps),] )
  num.point = d                     # suppose we have 'num.point' particles
  thABC = rnorm(num.point, m, sqrt(S))
  
  # Select a smaller rectangular domain D_2 in D_1
  Domain = range(thABC)
  th = matrix(0,num.point,p)
  
  # Generate d number of particles over the domain
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2]-Domain[1])*th[,1]+Domain[1]
  
  # Generate auxiliary variables at the sample mean of the particles
  hat = mean(th)
  Sample = pAuxStatGibb(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:num.point){
    cal = Sample*(th[i]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  if(LikEm){
    lhX = c()
    for(i in 1:num.point){ lhX[i] = stat*th[i] }
    y = lhX - y
  }
  
  # Fit a Gaussian process
  m = km(~ ., design = matrix(th[1:num.point], ncol=1), response = matrix(y[1:num.point], ncol=1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = hat
  cov = var(thABC)
  
  # calculating an initial value for the normalizing ftn / full likelihood
  l.h.X = hat*stat
  
  # Kriging
  x.point = data.frame("design" = theta[1])
  # x.point = data.frame("design" = 0.43)
  pred.m = predict(m, x.point, "UK")
  lhXZ = ypred = pred.m$mean
  
  
  if(LikEm){
    result = GPmcmcLik(Niter, theta, sqrt(cov), lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  } else {
    result = GPmcmcNorm(Niter, theta, sqrt(cov), lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  }
  
  return(result)
}


### Gaussian Process MCMC with ABC particle generation
# IS approximation with independent samples
IsingGPmcmcPrior = function(Niter, d, N, X, stat, cycle = 100, LikEm = FALSE, num){ 
  p = 1
  Domain = c(0, 1)
  num.point = 3000           
  th = matrix(0,num.point,p)
  
  # Generate D Latin hypercube design points over the domain D_1
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2]-Domain[1])*th[,1]+Domain[1]
  
  # Simulate the auxiliary variable for each design point
  summat = pAuxStatGibb2(X, 100, th, num) # 100 cycle
  
  # Choose some of them based on the distance b/w simulated and observed summary statistics
  dist = sqrt( ( rep(stat, num.point) - summat[1:num.point,] )^2 )
  eps = quantile(dist, probs = 0.03)
  m =  mean(th[which(dist < eps),])
  S = var( th[which(dist < eps),] )
  num.point = d                     # suppose we have 'num.point' particles
  thABC = rnorm(num.point, m, sqrt(S))
  
  # Select a smaller rectangular domain D_2 in D_1
  Domain = range(thABC)
  th = matrix(0,num.point,p)
  
  # Generate d number of particles over the domain
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2]-Domain[1])*th[,1]+Domain[1]
  
  # Generate auxiliary variables at the sample mean of the particles
  hat = mean(th)
  Sample = pAuxStatGibb(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:num.point){
    cal = Sample*(th[i]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  if(LikEm){
    lhX = c()
    for(i in 1:num.point){ lhX[i] = stat*th[i] }
    y = lhX - y
  }
  
  # Fit a Gaussian process
  m = km(~ ., design = matrix(th[1:num.point], ncol=1), response = matrix(y[1:num.point], ncol=1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = hat
  cov = var(thABC)
  
  # calculating an initial value for the normalizing ftn / full likelihood
  l.h.X = hat*stat
  
  # Kriging
  x.point = data.frame("design" = theta[1])
  # x.point = data.frame("design" = 0.43)
  pred.m = predict(m, x.point, "UK")
  lhXZ = ypred = pred.m$mean
  
  
  if(LikEm){
    result = GPmcmcLikPrior(Niter, theta, sqrt(cov), lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  } else {
    result = GPmcmcNormPrior(Niter, theta, sqrt(cov), lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  }
  
  return(result)
}



# ------------------------------------------------------------------------------
# Approximate intractable terms in our diagnostics
# ------------------------------------------------------------------------------
### approximate score function
f_uhat_b = function(Sx, Sy){
  return(Sx - mean(Sy))
}

### approximate H matrix
f_dhat_b = function(Sx, Sy){
  Sybar = mean(Sy)
  uhat = Sx - Sybar
  return( -mean(Sy^2) + Sybar * Sybar + uhat * uhat)
}

f_Vhat_bm = function(D, N){
  n = length(D)
  b = floor(min(n^(1/3), N^(1/3)))
  a = floor(n/b)

  dbarbk = sapply(1:a, function(k) return(mean(D[((k - 1) * b + 1):(k * b)])))
  dbar = mean(dbarbk)
  sigmahat = b * sum((dbarbk - dbar)^2) / (a-1)

  return( sigmahat )
}

ACD_bm = function(D, N){
  n = length(D)
  dbar = mean(D)
  
  res = n * (dbar^2) / f_Vhat_bm(D, N)
  return(res)
}



