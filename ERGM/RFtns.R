### Gaussian Process MCMC with ABC particle generation
ergmGPmcmc = function(Niter, d, N, X, stat, Domain, cycle = 100, LikEm = FALSE, num){ 
  p = 2
  # num.point = 3000           
  num.point = d*10 # LikEm2, LikEm3       
  th = matrix(0,num.point,p)
  
  # Generate D Latin hypercube design points over the domain D_1
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] = (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  
  # Simulate the auxiliary variable for each design point
  summat = pAuxSamp(X, cycle, th, 1, num) # 100 cycle
  
  # Choose some of them based on the distance b/w simulated and observed summary statistics
  dist = sqrt(  apply(  ( matrix( rep(stat, num.point), num.point, byrow = T) - summat[1:num.point,,1] )^2, 1, sum ) )
  eps = quantile(dist, probs = 0.03)
  m =  apply(  th[which(dist < eps),], 2, mean)
  S = cov( th[which(dist < eps),] )
  num.point = d                     # suppose we have 'num.point' particles
  thABC = mvrnorm(num.point, m, S)
  
  # Select a smaller rectangular domain D_2 in D_1
  Domain = rbind(apply(thABC,2,min), apply(thABC,2,max))  
  th = matrix(0,num.point,p)
  
  # Generate d number of particles over the domain
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] = (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  
  # Generate auxiliary variables at the sample mean of the particles
  hat = apply(thABC, 2, mean)
  Sample = pResponseErgm(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:num.point){
    cal = Sample%*%(th[i,]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  
  if(LikEm){
    lhX = c()
    for(i in 1:num.point){ lhX[i] = stat%*%th[i,] }
    y = lhX - y
  } 
  
  # Fit a Gaussian process
  m = km(~ ., design = th[1:num.point,], response = matrix(y[1:num.point],num.point,1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = matrix(hat,1)
  COV = cov(thABC)
  
  # calculating an initial value for the normalizing ftn / full likelihood
  # Kriging
  x.point = data.frame( X1 = theta[1,1], X2 = theta[1,2] )
  pred.m = predict(m, x.point, "UK")
  lhXZ = pred.m$mean
  
  if(LikEm){
    result = GPmcmcErgmLik(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  } else {
    result = GPmcmcErgmNorm(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  }
  
  return(result)
}


### Gaussian Process MCMC with DMH particle generation
ergmGPmcmcDMH = function(Niter, d, N, X, stat, parameter, hat, cycle = 100, LikEm = FALSE, num){ 
  th = parameter[1:d,]
  
  # Generate auxiliary variables at the sample mean of the particles
  Sample = pResponseErgm(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:d){
    cal = Sample%*%(th[i,]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  
  if(LikEm){
    lhX = c()
    for(i in 1:d){ lhX[i] = stat%*%th[i,] }
    y = lhX - y
  } 
  
  # Fit a Gaussian process
  m = km(~ ., design = th[1:d,], response = matrix(y[1:d],d,1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = matrix(hat,1)
  COV = diag(c(0.1, 0.1))
  
  # calculating an initial value for the normalizing ftn / full likelihood
  # Kriging
  x.point = data.frame( X1 = theta[1,1], X2 = theta[1,2] )
  pred.m = predict(m, x.point, "UK")
  lhXZ = pred.m$mean
  
  if(LikEm){
    result = GPmcmcErgmLik(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  } else {
    result = GPmcmcErgmNorm(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  }
  
  return(result)
}


### Gaussian Process MCMC with ABC particle generation
ergmABCLikPrior = function(Niter, d, N, X, stat, Domain, cycle = 100, num){ 
  p = 2
  # num.point = 3000           
  num.point = d*10 # LikEm2, LikEm3       
  th = matrix(0,num.point,p)
  
  # Generate D Latin hypercube design points over the domain D_1
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] = (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  
  # Simulate the auxiliary variable for each design point
  summat = pAuxSamp(X, cycle, th, 1, num) # 100 cycle
  
  # Choose some of them based on the distance b/w simulated and observed summary statistics
  dist = sqrt(  apply(  ( matrix( rep(stat, num.point), num.point, byrow = T) - summat[1:num.point,,1] )^2, 1, sum ) )
  eps = quantile(dist, probs = 0.03)
  m =  apply(  th[which(dist < eps),], 2, mean)
  S = cov( th[which(dist < eps),] )
  num.point = d                     # suppose we have 'num.point' particles
  thABC = mvrnorm(num.point, m, S)
  
  # Select a smaller rectangular domain D_2 in D_1
  Domain = rbind(apply(thABC,2,min), apply(thABC,2,max))  
  th = matrix(0,num.point,p)
  
  # Generate d number of particles over the domain
  Design = lhsDesign(n=num.point, dimension=p, randomized=TRUE, seed=1)
  th = Design$design
  th[,1] = (Domain[2,1]-Domain[1,1])*th[,1]+Domain[1,1]
  th[,2] = (Domain[2,2]-Domain[1,2])*th[,2]+Domain[1,2]
  
  # Generate auxiliary variables at the sample mean of the particles
  hat = apply(thABC, 2, mean)
  Sample = pResponseErgm(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:num.point){
    cal = Sample%*%(th[i,]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  
  lhX = c()
  for(i in 1:num.point){ lhX[i] = stat%*%th[i,] }
  y = lhX - y
  
  # Fit a Gaussian process
  m = km(~ ., design = th[1:num.point,], response = matrix(y[1:num.point],num.point,1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = matrix(hat,1)
  COV = cov(thABC)
  
  # calculating an initial value for the normalizing ftn / full likelihood
  # Kriging
  x.point = data.frame( X1 = theta[1,1], X2 = theta[1,2] )
  pred.m = predict(m, x.point, "UK")
  lhXZ = pred.m$mean
  
  result = GPmcmcErgmLikPrior(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  
  return(result)
}


### Gaussian Process MCMC with DMH particle generation
ergmDMHLikPrior = function(Niter, d, N, X, stat, parameter, hat, cycle = 100, num){ 
  th = parameter[1:d,]
  
  # Generate auxiliary variables at the sample mean of the particles
  Sample = pResponseErgm(X, cycle, hat, N, num) # summary statistics
  
  # Calculate importance sampling approximation to log normalizing function
  y = c()
  for(i in 1:d){
    cal = Sample%*%(th[i,]-hat)
    mx = max(cal)
    y[i] = mx + log( mean( exp(cal-mx) ) )
  }
  
  lhX = c()
  for(i in 1:d){ lhX[i] = stat%*%th[i,] }
  y = lhX - y 
  
  # Fit a Gaussian process
  m = km(~ ., design = th[1:d,], response = matrix(y[1:d],d,1), covtype = "matern3_2")
  
  # GLS beta coefficients
  beta.hat.gls = c(coef(m)$trend1,coef(m)$trend2,coef(m)$trend3)
  range.hat = coef(m)$range
  sd2.hat = coef(m)$sd2
  
  # initial theta and covariance matrix for proposal
  theta = matrix(hat,1)
  COV = diag(c(0.1, 0.1))
  
  # calculating an initial value for the normalizing ftn / full likelihood
  # Kriging
  x.point = data.frame( X1 = theta[1,1], X2 = theta[1,2] )
  pred.m = predict(m, x.point, "UK")
  lhXZ = pred.m$mean
  
  result = GPmcmcErgmLikPrior(Niter, theta, COV, lhXZ, beta.hat.gls, c(range.hat, sd2.hat), th, y, stat)[-1,]
  
  return(result)
}


### approximate score function
f_Uhat = function(Sx, Sy){
  return( as.vector(Sx - colMeans(Sy)) )
}


### approximate H matrix
f_Hhat = function(Sx, Sy){
  
  eS1 = mean(Sy[,1])
  eS2 = mean(Sy[,2])
  eS1S1 = mean(Sy[,1] * Sy[,1])
  eS2S2 = mean(Sy[,2] * Sy[,2])
  eS1S2 = mean(Sy[,1] * Sy[,2])
  
  H11 = - eS1S1 + eS1 * eS1
  H12 = - eS1S2 + eS1 * eS2
  H22 = - eS2S2 + eS2 * eS2
  
  return( as.vector(c(H11, H12, H22)) )
}


### approximate J matrix
f_Jhat = function(Sx, Sy){
  
  eS1 = mean(Sy[,1])
  eS2 = mean(Sy[,2])
  
  J11 = Sx[1] * Sx[1] -  2 * Sx[1] * eS1 + eS1 * eS1
  J12 = Sx[1] * Sx[2] - Sx[1] * eS2 - Sx[2] * eS1 + eS1 * eS2
  J22 = Sx[2] * Sx[2] - 2 * Sx[2] * eS2 + eS2 * eS2
  
  return( as.vector(c(J11, J12, J22)) )
}


### approximate A matrix
f_Ahat = function(Sx, Sy){
  
  eS1 = mean(Sy[,1])
  eS2 = mean(Sy[,2])
  eS1S1 = mean(Sy[,1] * Sy[,1])
  eS2S2 = mean(Sy[,2] * Sy[,2])
  eS1S2 = mean(Sy[,1] * Sy[,2])
  
  eS1S1S1 = mean(Sy[,1] * Sy[,1] * Sy[,1])
  eS1S1S2 = mean(Sy[,1] * Sy[,1] * Sy[,2])
  eS1S2S2 = mean(Sy[,1] * Sy[,2] * Sy[,2])
  eS2S2S2 = mean(Sy[,2] * Sy[,2] * Sy[,2])
  
  A11 = - 2 * Sx[1] * ( eS1S1 -  eS1 * eS1 ) + 5 * eS1 * eS1S1 - 4 * eS1 * eS1 * eS1 - eS1S1S1
  A12 = - 2 * Sx[1] * ( eS1S2 - eS1 * eS2 ) + 4 * eS1 * eS1S2 + eS2 * eS1S1 - 4 * eS1 * eS1 * eS2 - eS1S1S2
  A21 = - Sx[2] * ( eS1S1 - eS1 * eS1 ) - Sx[1] * ( eS1S2 - eS1 * eS2 ) + 3 * eS1 * eS1S2 + 2 * eS2 * eS1S1 - 4 * eS1 * eS1 * eS2 - eS1S1S2
  A22 = - Sx[1] * ( eS2S2 - eS2 * eS2 ) - Sx[2] * ( eS1S2 - eS1 * eS2 ) + 3 * eS2 * eS1S2 + 2 * eS1 * eS2S2 - 4 * eS1 * eS2 * eS2 - eS1S2S2
  A31 = - 2 * Sx[2] * ( eS1S2 - eS1 * eS2 ) + 4 * eS2 * eS1S2 + eS1 * eS2S2 - 4 * eS1 * eS2 * eS2 - eS1S2S2
  A32 = - 2 * Sx[2] * ( eS2S2 -  eS2 * eS2 ) + 5 * eS2 * eS2S2 - 4 * eS2 * eS2 * eS2 - eS2S2S2
  
  return( as.vector(c(A11, A12, A21, A22, A31, A32)) )
}


### biased but consistent estimate of V
f_Vhat = function(U, H, J, A){
  n = nrow(U)
  d = J + H
  Abar = matrix(colMeans(A), 3, 2, byrow = T)
  Bbar = (t(d) %*% d) / n
  Cbar = (t(U) %*% d) / n
  Hbar = matrix(colMeans(cbind(H[,1], H[,2], H[,2], H[,3])), 2, 2)
  Jbar = matrix(colMeans(cbind(J[,1], J[,2], J[,2], J[,3])), 2, 2)
  invHhar = solve(Hbar)
  
  return( Bbar - 2 * Abar %*% invHhar %*% Cbar + Abar %*% invHhar %*% Jbar %*% invHhar %*% t(Abar) )
}


### ACD using the biased but consistent estimate of V
ACD = function(U, H, J, A){
  n = nrow(U)
  dbar = colMeans(J + H)
  
  res = n * t(dbar) %*% solve(f_Vhat(U, H, J, A)) %*% dbar
  return(as.vector(res))
}






