// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>

#define minimum(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// =============================================================================
//  sampling (rejection sampling) from COMP distribution
// =============================================================================
// Evaluate unnormalised density of the COM-poisson distribution 
// if fudge is set to lgammafn(mode+1) then the unnormalised density is one at the mode
// If mode and fudge are set to 0 then the usual unnormalised density is computed
// [[Rcpp::export]]
double unnorm_ldcpois(double x, double mu, double nu, double mode, double fudge) {
  return nu*((x-mode)*log(mu)-lgamma(x+1)+fudge);
}


// Sample from a geometric distribution truncated to {0,1,...,n}
// u is a U[0,1] realisation
// [[Rcpp::export]]
double truncated_geo_sample(double u, double logq, double n) {
  double C;
  if(logq > -DBL_EPSILON){
    return 0;
  } else {
    C = -expm1(logq*(n+1));
    return floor(log(1-C*u)/logq); 
  }
}


// Sample from a geometric distribution with range {0,1,2,...}
// u is a U[0,1] realisation
// [[Rcpp::export]]
double untruncated_geo_sample(double u, double logq) {
  if (logq > -DBL_EPSILON){
    return 0;
  } else {
    return floor(log(u)/logq); 
  }
}


// Compute the finite geometric series (1+q+...+q^deltax)
// [[Rcpp::export]]
double truncated_lweights(double deltax, double logq) {
  if (logq > -DBL_EPSILON)
    return log(deltax+1)+logq;
  return log1p(-exp((deltax+1)*logq)) - log1p(-exp(logq));
}


// Compute the geometric series (1+q+...)
// [[Rcpp::export]]
double untruncated_lweights(double logq) {
  return -log1p(-exp(logq));
}


// Sample from an COMP distribution
// [[Rcpp::export]]
vec rCOMP(int n, double mu, double nu){ // mu = mode
  double negInf = -std::numeric_limits<float>::infinity();;
  
  double logmu, lmode, rmode, fudge, sd, lsd, rsd, maxlweight, logprob, x, u;
  int i, attempts;
  vec ldens(4), lweights(4), logq(4), sweights(4), result(n);
  logmu = log(mu);
  
  // Figure out mode and standard deviation
  lmode = ceil(mu)-1;
  fudge = lgamma(lmode+1);
  rmode = lmode+1;
  fudge = lgamma(lmode+1);
  sd = ceil(sqrt(mu)/sqrt(nu));
  if (sd < 5) {
    sd = 5;
  }
  
  // Set up two points at mode +/- sd
  lsd = round(lmode-sd);
  if (lsd < 0){
    lsd = -1;
  }
  rsd = round(rmode+sd);
  
  // Left most tail
  if (lsd == -1) {
    lweights[0] = negInf;
    logq[0] = 0;
    ldens[0] = negInf;
  } else {
    ldens[0] = unnorm_ldcpois(lsd, mu, nu, lmode, fudge);
    if (lsd == 0) {
      lweights[0] = ldens[0];
      logq[0] = 0;
    } else {
      logq[0] = nu * (-logmu + log(lsd));
      lweights[0] = ldens[0] + truncated_lweights(lsd, logq[0]);
    }
  }
  
  // within 1sd to the left of the mode
  ldens[1] = 0;
  if (lmode == 0) {
    lweights[1] = 0;
    logq[1] = 1;
  } else {
    logq[1] = nu * (-logmu + log(lmode));
    lweights[1] = truncated_lweights(lmode-lsd-1, logq[1]);
  }
  
  // within 1sd to the right of the mode
  logq[2] = nu * (logmu - log(rmode+1));
  ldens[2] = nu * (logmu - log(rmode));
  lweights[2] = ldens[2] + truncated_lweights(rsd-rmode-1, logq[2]);
  
  // right tail
  logq[3] = nu * (logmu - log(rsd+1));
  ldens[3] = unnorm_ldcpois(rsd, mu, nu, lmode, fudge);
  lweights[3] = ldens[3] + untruncated_lweights(logq[3]);
  
  // Find maximum log-weight
  maxlweight = lweights[0];
  for (i = 1; i < 4; i++){
    if (lweights[i] > maxlweight) { maxlweight = lweights[i]; }
  }
  // Compute the cumulative sum of the weights
  for (i = 0; i < 4; i++) {
    lweights[i] = lweights[i]-maxlweight;
    sweights[0] = exp(lweights[0]);
  }
  for (i = 1; i < 4; i++) {
    sweights[i]=sweights[i-1]+exp(lweights[i]);
  }
  
  // Draw the sample by rejection sampling
  attempts = 0;
  for (i = 0; i < n; i++) {
    while (TRUE) {
      attempts = attempts + 1;
      u = randu() * sweights[3];
      if (u < sweights[0]) {
        u = u / sweights[0];
        x = truncated_geo_sample(u, logq[0], lsd);
        logprob = ldens[0]+x*logq[0];
        x = lsd-x;
      } else {
        if (u < sweights[1]) {
          u = (u-sweights[0])/(sweights[1]-sweights[0]);
          x = truncated_geo_sample(u, logq[1], lmode-lsd-1);
          logprob = ldens[1]+x*logq[1];
          x = lmode - x;
        } else {
          if (u<sweights[2]) {
            u = (u-sweights[1])/(sweights[2]-sweights[1]);
            x = truncated_geo_sample(u, logq[2], rsd-rmode-1);
            logprob = ldens[2]+x*logq[2];
            x = rmode + x;
          } else {
            u = (u-sweights[2])/(sweights[3]-sweights[2]);
            x = untruncated_geo_sample(u, logq[3]);
            logprob = ldens[3]+x*logq[3];
            x = rsd + x;
          }
        }
      }
      if (log(randu()) < unnorm_ldcpois(x, mu, nu, lmode, fudge) - logprob) {
        result[i] = x;
        break;
      }
    }
  }
  return result;
}



// Auxiliary variable generation in parallel
// [[Rcpp::export]]
mat rCOMP_parallel(int n, vec mode, double nu, int num){
  int N = mode.size();
  mat aux = zeros(n, N);
  
  #pragma omp parallel num_threads(num)
  {
  #pragma omp for
    for(int i = 0; i < N; i++){
      aux.col(i) = rCOMP(n, mode[i], nu);
    }
  }
  
  return aux;
}


// =============================================================================
// Distribution functions
// =============================================================================

// Unnormalized log likelihood of COMP
// [[Rcpp::export]]
double modeCOMP_logh(vec y, vec mode, double nu){
  double result = sum( nu * ( y % log(mode) -  lgamma(y+1) ) );
  return result;
}

// Unnormalized likelihood of COMP
// [[Rcpp::export]]
double modeCOMP_h(vec y, vec mode, double nu){
  double result = sum( nu * ( y % log(mode) -  lgamma(y+1) ) );
  return exp(result);
}


// normalizing function of COMP
// [[Rcpp::export]]
double modeCOMP_Z(double mode, double nu, int summax){
  vec y = regspace(0, summax);
  vec x = nu * ( y * log(mode) - lgamma(y + 1) );
  double result = x.max() + log(sum(exp(x - x.max())));
  return exp(result);
}


// normalizing function of COMP in parallel
// [[Rcpp::export]]
vec modeCOMP_Z_parallel(vec mode, double nu, int summax, int num){
  int N = mode.size();
  vec result(N);
#pragma omp parallel num_threads(num)
{
#pragma omp for
  for(int i = 0; i < N; i++){
    result[i] = modeCOMP_Z(mode[i], nu, summax);
  }
}
return result;
}



// log of normalizing function of COMP
// [[Rcpp::export]]
double modeCOMP_logZ(double mode, double nu, int summax){
  vec y = regspace(0, summax);
  vec x = nu * ( y * log(mode) - lgamma(y + 1) );
  double result = x.max() + log(sum(exp(x - x.max())));
  return result;
}


// log of normalizing function of COMP in parallel
// [[Rcpp::export]]
vec modeCOMP_logZ_parallel(vec mode, double nu, int summax, int num){
  int N = mode.size();
  vec result(N);
#pragma omp parallel num_threads(num)
{
#pragma omp for
  for(int i = 0; i < N; i++){
    result[i] = modeCOMP_logZ(mode[i], nu, summax);
  }
}
return result;
}


// log likelihood of COMP
// [[Rcpp::export]]
vec modeCOMP_logd(vec y, vec mode, double nu, int summax, int num){
  vec result = nu * ( y % log(mode) -  lgamma(y+1) ) - modeCOMP_logZ_parallel(mode, nu, summax, num);
  return result;
}


// Unnormalized log normal
// [[Rcpp::export]]
double Normal_logh_mar(double y, double mu, double sig2){
  double result = - 0.5 * pow(y - mu, 2) / sig2;
  return result;
} 

// Unnormalized joint log normal
// [[Rcpp::export]]
double Normal_logh(vec y, vec mu, vec sig2){
  double result = sum(- 0.5 * (y - mu) % (y - mu) / sig2);
  return result;
} 


// Unnormalized log multivariate normal
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
} 


// log Poisson density
// [[Rcpp::export]]
double Poisson_logd(double y, double lambda){
  double result = log(y) * log(lambda) - lambda - lgamma(y + 1);
  return result;
} 


// =============================================================================
// COMP Regression Models
// =============================================================================
// Approximation method for inference
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List appxCOMPreg(int outer, vec y, mat X, vec beta, double alpha,
                 List block, vec sigma2, List COV,
                 bool updateCOV, int updateUntil,
                 int adaptInterval, double adaptFactorExponent,
                 int summax, int thin, int num){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj, adapIter = 1;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1);
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    cholCOV[i] = trans( chol( sigma2[i] * as<arma::mat>(COV[i]) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = exp(alpha);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){

        if( (k >= adaptInterval) && (k <= updateUntil) && (k - (adaptInterval * trunc(k / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k-adaptInterval, k-1) ) / adaptInterval;
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(0, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if(j == m1-1){ adapIter = adapIter + 1; }
        }
      } 
      
      
      vec dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = exp(alphaprop);
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        sum(modeCOMP_logZ_parallel(mode, nu, summax, num)) - sum(modeCOMP_logZ_parallel(modeprop, nuprop, summax, num)) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      }
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    }
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV);
}


// known dispersion
// [[Rcpp::export]]
List appxCOMPreg_knu(int outer, vec y, mat X, vec beta, double nu,
                     List block, vec sigma2, List COV, bool updateCOV,
                     int adaptInterval, double adaptFactorExponent, int adapIter,
                     int summax, int thin, int num){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p), posterior_thined, accprob = zeros(outer, m1), postj;
  vec betaprop(p), mode(N), modeprop(N), par(p), parprop(p);
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  par.elem(index_beta) = beta;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        }  
      } 
      
      
      vec dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      } 
      
      betaprop = parprop.elem(index_beta);
      
      modeprop = exp( X * betaprop );
      
      logprob = modeCOMP_logh(y, modeprop, nu) - modeCOMP_logh(y, mode, nu) +
        sum(modeCOMP_logZ_parallel(mode, nu, summax, num)) - sum(modeCOMP_logZ_parallel(modeprop, nu, summax, num)) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        mode = modeprop;
        
        accprob(k, j) = 1;
      } 
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    } 
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    } 
    
  } 
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
} 


// -----------------------------------------------------------------------------
// Exchange algorithm for inference
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List exactCOMPreg(int outer, vec y, mat X, vec beta, double alpha,
                 List block, vec sigma2, List COV, bool updateCOV,
                 int adaptInterval, double adaptFactorExponent, int adapIter,
                 int thin){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1), aux(N), dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = exp(alpha);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        }
      } 
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = exp(alphaprop);
      
      // aux = trans( rCOMP_parallel(1, modeprop, nuprop, num) );
      for(int i = 0; i < N; i ++){
        aux[i] = ( rCOMP(1, modeprop[i], nuprop) )[0];
      }
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nuprop) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      }
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    }
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
}


// known dispersion
// [[Rcpp::export]]
List exactCOMPreg_knu(int outer, vec y, mat X, vec beta, double nu,
                      List block, vec sigma2, List COV, bool updateCOV,
                      int adaptInterval, double adaptFactorExponent, int adapIter,
                      int thin){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p), parprop(p), aux(N), dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  par.elem(index_beta) = beta;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        } 
      } 
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      } 
      
      betaprop = parprop.elem(index_beta);
      
      modeprop = exp( X * betaprop );
      
      // aux = trans( rCOMP_parallel(1, modeprop, nu, num) );
      for(int i = 0; i < N; i ++){
        aux[i] = ( rCOMP(1, modeprop[i], nu) )[0];
      }
      
      logprob = modeCOMP_logh(y, modeprop, nu) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nu) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        mode = modeprop;
        
        accprob(k, j) = 1;
      }
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    } 
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    } 
    
  } 
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
} 


// -----------------------------------------------------------------------------
// Exchange with auxiliary sampling via truncation
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
int sample_cpp(vec x, vec probs){
  double u = randu();
  int nn = probs.size();
  vec csprobs = zeros(nn+1);
  csprobs.rows(1,nn) = cumsum(probs);
  csprobs.insert_rows(0, 0);
  int i = 0, res;
  while(u > csprobs[i]){
    res = x[i];
    i = i + 1;
  }
  return res;
}  


// [[Rcpp::export]]
int rcomp_trunc(double mode, double nu, int summax){
  vec y = regspace(0, summax);
  vec x = nu * ( y * log(mode) - lgamma(y + 1) );
  double logZ = x.max() + log(sum(exp(x - x.max())));
  vec probs = exp(nu * ( y * log(mode) - lgamma(y+1) ) - logZ);
  int auxy = sample_cpp(y, probs);
  return auxy;
} 


// [[Rcpp::export]]
vec rcomp_trunc_parallel(vec mode, double nu, int summax, int num){
  int N = mode.size();
  vec y = regspace(0, summax), res(N);
  
#pragma omp parallel num_threads(num)
{
#pragma omp for
  for(int i = 0; i < N; i++){
    // vec y = regspace(0, summax);
    vec x = nu * ( y * log(mode[i]) - lgamma(y + 1) );
    double logZ = x.max() + log(sum(exp(x - x.max())));
    vec probs = exp(nu * ( y * log(mode[i]) - lgamma(y+1) ) - logZ);
    res[i] = sample_cpp(y, probs);
  }
} 

return res; 
}



// [[Rcpp::export]]
List compExchangeTrunc(int outer, vec y, mat X, vec beta, double alpha,
                       List block, vec sigma2, List COV, bool updateCOV,
                       int adaptInterval, double adaptFactorExponent, int adapIter,
                       int summax, int thin, int num){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1), aux, dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }  
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = exp(alpha);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        } 
      }  
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }  
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = exp(alphaprop);
      
      aux = rcomp_trunc_parallel(modeprop, nuprop, summax, num);
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nuprop) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      } 
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      }
    } 
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    } 
    
  } 
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
} 



// known dispersion
// [[Rcpp::export]]
List compExchangeTrunc_knu(int outer, vec y, mat X, vec beta, double nu,
                           List block, vec sigma2, List COV, bool updateCOV,
                           int adaptInterval, double adaptFactorExponent, int adapIter,
                           int summax, int thin, int num){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj;
  mat posterior(outer, p), posterior_thined, accprob = zeros(outer, m1), prej, postj;
  vec betaprop(p), mode(N), modeprop(N), par(p), parprop(p), aux, dummy;
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  par.elem(index_beta) = beta;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    uvec index = block[i];
    mj = index.size();
    cholCOV[i] = trans( chol( sigma2[i] * ( as<arma::mat>(COV[i]) + 0.001 * diagmat(ones(mj)) ) ) );
  }
  
  
  // initial parameter values
  mode = exp(X * beta);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k+1 >= adaptInterval) && (k+1 - (adaptInterval * trunc((k+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k+1-adaptInterval, k-1) ) / (adaptInterval-1);
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(k+1-adaptInterval, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if( j == m1-1 ){ adapIter = adapIter + 1; }
        } 
      }  
      
      
      dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      }  
      
      betaprop = parprop.elem(index_beta);
      
      modeprop = exp( X * betaprop );
      
      aux = rcomp_trunc_parallel(modeprop, nu, summax, num);
      
      logprob = modeCOMP_logh(y, modeprop, nu) - modeCOMP_logh(y, mode, nu) +
        modeCOMP_logh(aux, mode, nu) - modeCOMP_logh(aux, modeprop, nu) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        mode = modeprop;
        
        accprob(k, j) = 1;
      } 
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      } 
    } 
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    }  
    
  }  
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COV") = COV);
}  


// -----------------------------------------------------------------------------
// Noisy exchange algorithm for inference
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
vec ratio_modeCOMP_h_parallel(mat y, vec mode, double nu, vec modeprop, double nuprop, int num){
  int N = y.n_cols;
  vec res = zeros(N);
  
  #pragma omp parallel num_threads(num)
  {
  #pragma omp for
    for(int i = 0; i < N; i++){
      res[i] = modeCOMP_h(y.col(i), mode, nu) / modeCOMP_h(y.col(i), modeprop, nuprop);
    }
  }
  return res;
}


// [[Rcpp::export]]
List noisyCOMPreg(int outer, int Ne, vec y, mat X, vec beta, double alpha,
                  List block, vec sigma2, List COV,
                  bool updateCOV, int updateUntil,
                  int adaptInterval, double adaptFactorExponent,
                  int thin, int num){
  
  int p = beta.size(), N = y.size(), iter = 0, m1 = block.size(), mj, adapIter = 1;
  mat posterior(outer, p+1), posterior_thined, accprob = zeros(outer, m1), prej, postj, aux;
  vec betaprop(p), mode(N), modeprop(N), par(p+1), parprop(p+1);
  vec rhat = zeros(m1), gamma1 = zeros(m1), gamma2 = zeros(m1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234;
  double logprob, u, alphaprop, nu, nuprop;
  uvec index_beta(p);
  for(int i = 0; i < p; i++){ index_beta[i] = i; }
  int index_alpha = p;
  par.elem(index_beta) = beta;
  par[index_alpha] = alpha;
  
  List cholCOV(m1);
  for(int i = 0; i < m1; i++){
    cholCOV[i] = trans( chol( sigma2[i] * as<arma::mat>(COV[i]) ) );
  } 
  
  
  // initial parameter values
  mode = exp(X * beta);
  nu = exp(alpha);
  
  // Start of MCMC Chain
  for(int k = 0; k < outer; k++){
    
    for(int j = 0; j < m1; j ++){
      
      uvec index = block[j];
      mj = index.size();
      
      // update proposal distribution
      if( updateCOV ){
        
        if( (k >= adaptInterval) && (k <= updateUntil) && (k - (adaptInterval * trunc(k / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(j);
          rhat[j] = sum( dummyaccprob.rows(k-adaptInterval, k-1) ) / adaptInterval;
          gamma1[j] = 1 / pow(adapIter, c1);
          gamma2[j] = c0 * gamma1[j];
          sigma2[j] = exp( log(sigma2[j]) + gamma2[j] * (rhat[j] - ropt) );
          
          postj = posterior.cols(index);
          COV[j] = as<arma::mat>(COV[j]) + gamma1[j] * ( cov( postj.rows(0, k-1) ) - as<arma::mat>(COV[j]) );
          
          cholCOV[j] = trans( chol( sigma2[j] * ( as<arma::mat>(COV[j]) + 0.001 * diagmat(ones(mj)) ) ) );
          
          if(j == m1-1){ adapIter = adapIter + 1; }
        } 
      } 
      
      
      vec dummy = as<arma::mat>(cholCOV[j]) * randn(mj);
      parprop = par;
      for(int l = 0; l < mj; l ++){
        parprop[index[l]] = par[index[l]] + dummy[l];
      } 
      
      betaprop = parprop.elem(index_beta);
      alphaprop = parprop[index_alpha];
      
      modeprop = exp( X * betaprop );
      nuprop = exp(alphaprop);
      
      aux = trans( rCOMP_parallel(Ne, modeprop, nuprop, num) );
      
      logprob = modeCOMP_logh(y, modeprop, nuprop) - modeCOMP_logh(y, mode, nu) +
        log( mean( ratio_modeCOMP_h_parallel(aux, mode, nu, modeprop, nuprop, num) ) ) +
        MVN_logh(parprop.elem(index), zeros(mj), diagmat(ones(mj))/100) -
        MVN_logh(par.elem(index), zeros(mj), diagmat(ones(mj))/100);
      
      u = log( randu() );
      if( u < logprob ){
        par = parprop;
        beta = betaprop;
        alpha = alphaprop;
        mode = modeprop;
        nu = nuprop;
        
        accprob(k, j) = 1;
      } 
      
      for(int l = 0; l < mj; l ++){
        posterior(k,index[l]) = par[index[l]];
      } 
    }
    
    // thining
    if(k - (thin * trunc(k / thin)) == 0){
      iter = iter + 1;
      posterior_thined.insert_rows(iter-1, posterior.row(k));
    } 
    
  } 
  
  return Rcpp::List::create(Rcpp::Named("Sample") = posterior_thined,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("COV") = COV);
} 


// =============================================================================
// ACD
// =============================================================================

// [[Rcpp::export]]
mat rCOMP2_parallel(mat X, vec theta, int N, int num){
  int n = X.n_rows, p = theta.size();
  vec mode, beta = theta.rows(0, p-2);
  mat aux(N, n), sumstat(N, p);
  double alpha = theta[p-1], nu;

  mode = exp(X * beta);
  nu = exp(alpha);

  aux = rCOMP_parallel(N, mode, nu, num);
  sumstat.cols(0, p-2) = aux * X;
  sumstat.col(p-1) = sum(lgamma(aux+1), 1); 
  
  return sumstat;
}


// [[Rcpp::export]]
vec rcppf_Uhat(vec Sx, mat Sy, vec theta){
  int p = theta.size(); 
  vec beta = theta.rows(0, p-2), eSy = mean(trans(Sy), 1), Uhat(p);
  vec Sx1 = Sx.rows(0, p-2), eSy1 = eSy.rows(0, p-2);
  double nu = exp(theta[p-1]);
  
  Uhat.rows(0, p-2) = nu * (Sx1 - eSy1);
  Uhat[p-1] = sum(beta % Sx1) - Sx[p-1] - sum(beta % eSy1) + eSy[p-1];
  
  return Uhat;
}


// [[Rcpp::export]]
mat rcppf_Hhat(vec Sx, mat Sy, vec theta){
  int p = theta.size(); 
  vec beta = theta.rows(0, p-2);
  mat Hhat(p, p);
  double nu = exp(theta[p-1]);
  
  for(int j = 0; j < p; j ++){
    for(int i = 0; i < p; i ++){
      
      if( (i < p-1) & (j < p-1) & (i <= j) ){
        
        Hhat(i,j) = Hhat(j,i) = - mean( nu * nu * Sy.col(i) % Sy.col(j) ) + mean(nu * Sy.col(i)) * mean(nu * Sy.col(j));
        
      } else if( (i < p-1) & (j == p-1) ) {
        
        Hhat(i,j) = Sx[i] - mean(Sy.col(i)) - 
          mean( nu * Sy.col(i) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) +
          mean( nu * Sy.col(i)) * mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) );
        
      } else if( (i == p-1) & (j < p-1) ) {
        
        Hhat(i,j) = Sx[j] - mean(Sy.col(j)) - 
          mean( nu * Sy.col(j) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) +
          mean( nu * Sy.col(j)) * mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) );
        
      } else if( (i == p-1) & (j == p-1) ) {
        
        Hhat(i,j) = - mean( ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) +
          mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) * mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) );
        
      }
    }  
  }
  
  return Hhat;
}



// // [[Rcpp::export]]
// vec rcppf_Ahat(vec Sx, mat Sy, vec theta, vec Uhat, mat Hhat){
//   int p = theta.size(), iter = 0; 
//   vec beta = theta.rows(0, p-2), Ahat(p * p * (p+1)/2);
//   double nu = exp(theta[p-1]);
//   
//   for(int m = 0; m < p; m ++){
//     for(int l = 0; l < p; l ++){
//       for(int k = 0; k <= l; k ++){
//         
//         if( (k < p-1) & (l < p-1) & (m < p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) + 
//             mean(nu * Sy.col(k)) * mean(nu * nu * Sy.col(l) % Sy.col(m)) +
//             mean(nu * Sy.col(l)) * mean(nu * nu * Sy.col(k) % Sy.col(m)) +
//             mean(nu * Sy.col(m)) * mean(nu * nu * Sy.col(k) % Sy.col(l)) -
//             mean(nu * nu * nu * Sy.col(k) % Sy.col(l) % Sy.col(m)) -
//             2 * mean(nu * Sy.col(k)) * mean(nu * Sy.col(l)) * mean(nu * Sy.col(m));
//           iter = iter + 1;
//           
//         } else if( (k < p-1) & (l < p-1) & (m == p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) +
//             mean(nu * Sy.col(k)) * ( mean(Sy.col(l)) + mean( nu * Sy.col(l) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) ) +
//             mean(nu * Sy.col(l)) * ( mean(Sy.col(k)) + mean( nu * Sy.col(k) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) ) +
//             mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) * mean( nu * nu *  Sy.col(k) % Sy.col(l) ) -
//             2 * mean( nu * Sy.col(k) % Sy.col(l) ) - 
//             mean( nu * nu * Sy.col(k) % Sy.col(l) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) -
//             2 * mean(nu * Sy.col(k)) * mean(nu * Sy.col(l)) * mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) );
//           iter = iter + 1;
//           
//         } else if( (k < p-1) & (l == p-1) & (m < p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) + 
//             mean(nu * Sy.col(k)) * ( mean(Sy.col(m)) + mean( nu * Sy.col(m) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) ) +
//             mean(nu * Sy.col(m)) * ( mean(Sy.col(k)) + mean( nu * Sy.col(k) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) ) +
//             mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) * mean( nu * nu * Sy.col(k) % Sy.col(m) ) -
//             2 * mean( nu * Sy.col(k) % Sy.col(m) ) - 
//             mean( nu * nu * Sy.col(k) % Sy.col(m) % ( Sy.cols(0, p-2) * beta - Sy.col(p-1) ) ) -
//             2 * mean(nu * Sy.col(k)) * mean(nu * Sy.col(m)) * mean( Sy.cols(0, p-2) * beta - Sy.col(p-1) );
//           iter = iter + 1;
//           
//         } else if( (k < p-1) & (l == p-1) & (m == p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) +
//             mean(nu * Sy.col(k)) * mean( (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) +
//             2 * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * ( mean(Sy.col(k)) + mean( nu * Sy.col(k) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) ) -
//             2 * mean( Sy.col(k) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             mean(nu * Sy.col(k) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             2 * mean(nu * Sy.col(k)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1));
//           iter = iter + 1;
//           
//         } else if( (k == p-1) & (l == p-1) & (m < p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) +
//             mean(nu * Sy.col(m)) * mean( (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) +
//             2 * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * ( mean(Sy.col(m)) + mean( nu * Sy.col(m) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) ) -
//             2 * mean( Sy.col(m) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             mean(nu * Sy.col(m) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             2 * mean(nu * Sy.col(m)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1));
//           iter = iter + 1;
//           
//         } else if( (k == p-1) & (l == p-1) & (m == p-1) ){
//           
//           Ahat[iter] = Uhat[k] * Hhat(l,m) + Uhat[l] * Hhat(k,m) +
//             3 * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * mean( (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             mean( (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) % (Sy.cols(0, p-2) * beta - Sy.col(p-1)) ) -
//             2 * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1)) * mean(Sy.cols(0, p-2) * beta - Sy.col(p-1));
//           iter = iter + 1;
//           
//         }
//       }
//     }
//   }
//   
//   return Ahat;
// }


// known dispersion ---------------------------------------
// [[Rcpp::export]]
mat rCOMP2_parallel_knu(mat X, vec beta, double nu, int N, int num){
  int n = X.n_rows, p = beta.size();
  vec mode;
  mat aux(N, n), sumstat(N, p);
  mode = exp(X * beta);

  aux = rCOMP_parallel(N, mode, nu, num);
  sumstat.cols(0, p-1) = aux * X;
  
  return sumstat;
}

// [[Rcpp::export]]
vec rcppf_Uhat_knu(vec Sx, mat Sy, vec beta, double nu){
  vec eSy = mean(trans(Sy), 1);
  vec Uhat = nu * (Sx - eSy);
  return Uhat;
}

// [[Rcpp::export]]
mat rcppf_Hhat_knu(vec Sx, mat Sy, double nu){
  int p = Sx.size();
  mat Hhat(p, p);

  for(int j = 0; j < p; j ++){
    for(int i = 0; i < p; i ++){

      if( i < j ){
        Hhat(i,j) = Hhat(j,i) = - mean( nu * nu * Sy.col(i) % Sy.col(j) ) + mean(nu * Sy.col(i)) * mean(nu * Sy.col(j));
      } else if( i == j ){
        Hhat(i,j) = - mean( nu * nu * Sy.col(i) % Sy.col(j) ) + mean(nu * Sy.col(i)) * mean(nu * Sy.col(j));
      } 
    }
  }

  return Hhat;
}
// -------------------------------------------------------



// [[Rcpp::export]]
vec AIKS(mat th, mat score, vec weight, double c, double beta, int i, int k, int num){
  int p = score.n_cols;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip;
  mat knot = zeros(k,p);
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    }
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    } 
  }
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    double wsq = pow(weight[i],2)*knot(i,j);
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        wsq += weight[i]*knot(m,j)*weight[m]*2;
      } 
    }
    H[j] = wsq;
  } 
  
  return(H);
} 





// [[Rcpp::export]]
// Compute w^2 for each row
mat pAIKS(mat th, mat score, vec weight, double c, double beta, int k, int num){
  int p = score.n_cols;
  mat H(k,p);
  omp_set_num_threads(num);
  
  int i;
#pragma omp parallel shared(H) private(i)
{ 
#pragma omp for schedule(static)
  for(i = 0; i < k; i++){
    vec res = AIKS(th, score, weight, c, beta, i, k, num);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  } 
}
return(H);
}



// [[Rcpp::export]]
vec AIKS_woWeight(mat th, mat score, double c, double beta, int i, int k, int num){
  int p = score.n_cols;
  double temp, temp2, bthi, bthip, k0, k0thi, k0thip, k0thithip;
  mat knot = zeros(k,p);
  
  for(int ip = i; ip < k; ip++){
    temp = pow(c, 2);
    for(int j = 0; j < p; j++){
      temp += pow(th(i,j) - th(ip,j), 2);
    } 
    
    for(int j = 0; j < p; j++){
      bthi = score(i,j);
      bthip = score(ip,j);
      temp2 = th(i,j) - th(ip,j);
      
      k0 = pow(temp, beta);
      k0thi = 2*beta*pow(temp, beta-1)*temp2;
      k0thip = -k0thi;
      k0thithip = -2*beta*pow(temp,beta-2)*(2*pow(temp2,2)*(beta-1)+temp);
      
      knot(ip,j) = bthi*bthip*k0 + bthi*k0thip + bthip*k0thi + k0thithip;
    }
  } 
  
  vec H = zeros(p);
  
  for(int j = 0; j < p; j++){
    double wsq = knot(i,j);
    if(i < k-1){
      for(int m = (i+1); m < k; m++){
        wsq += knot(m,j)*2;
      }
    } 
    H[j] = wsq;
  } 
  
  return(H);
} 





// [[Rcpp::export]]
// Compute w^2 for each row
mat pAIKS_woWeight(mat th, mat score, double c, double beta, int k, int num){
  int p = score.n_cols;
  mat H(k,p);
  omp_set_num_threads(num);
  
  int i;
#pragma omp parallel shared(H) private(i)
{ 
#pragma omp for schedule(static)
  for(i = 0; i < k; i++){
    vec res = AIKS_woWeight(th, score, c, beta, i, k, num);
    
    for(int j = 0; j < p; j++){
      H(i,j) = res[j];
    }
  }
}
return(H);
}





