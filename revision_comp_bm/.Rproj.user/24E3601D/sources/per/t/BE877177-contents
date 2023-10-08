f_Uhat_knu = function(Sx, Sy, nu){
  
  eSy = colMeans(Sy)
  return( nu * (Sx - eSy) )
  
}

f_dhat_knu = function(Sx, Sy, nu, Uhat){
  
  Hhat = rcppf_Hhat_knu(Sx, Sy, nu)
  Jhat = Uhat %*% t(Uhat)
  dhat = Hhat + Jhat
  
  return( dhat[upper.tri(dhat, diag = T)] )
}

f_Vhat_bm = function(d, N){
  n = nrow(d)
  b = floor(min(n^(1/3), N^(1/3)))
  a = floor(n/b)
  
  dbarbk = sapply(1:a, function(k) return(colMeans(d[((k - 1) * b + 1):(k * b),])))
  dbar = rowMeans(dbarbk)
  dummy = 0
  for(k in 1:a){
    dummy = dummy + (dbarbk[,k] - dbar) %*% t(dbarbk[,k] - dbar)
  }
  Sigmahat = b * dummy / (a-1)
  
  return( Sigmahat )
}


### ACD using the biased but consistent estimate of V
ACD_bm = function(d, N){
  n = nrow(d)
  dbar = colMeans(d)
  res = n * t(dbar) %*% solve(f_Vhat_bm(d, N)) %*% dbar
  return(as.vector(res))
}




