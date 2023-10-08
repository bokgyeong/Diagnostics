### Pseudolikelihood for MPLE
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


### Approximate score function
f_uhat_b = function(Sx, Sy){
  return(Sx - mean(Sy))
}


### Approximate difference between hessian and outerproject of score
f_dhat_b = function(Sx, Sy){
  Sybar = mean(Sy)
  uhat = Sx - Sybar
  return( -mean(Sy^2) + Sybar * Sybar + uhat * uhat)
}


### Approximate asymototic variance
f_Vhat_bm = function(D, N){
  n = length(D)
  b = floor(min(n^(1/3), N^(1/3)))
  a = floor(n/b)

  dbarbk = sapply(1:a, function(k) return(mean(D[((k - 1) * b + 1):(k * b)])))
  dbar = mean(dbarbk)
  sigmahat = b * sum((dbarbk - dbar)^2) / (a-1)

  return( sigmahat )
}

                  
### Compute ACD
ACD_bm = function(D, N){
  n = length(D)
  dbar = mean(D)
  
  res = n * (dbar^2) / f_Vhat_bm(D, N)
  return(res)
}
