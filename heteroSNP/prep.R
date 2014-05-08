library(VGAM)

lambda <- 0.5;

marginalLikelihood <- function(params, data) {
  
  alpha <- params[1];
  beta <- params[2];
  vec <-data;
  dlen <- length(vec);
  
  As <- vec[seq(1, dlen, 2)];
  Bs <- vec[seq(2, dlen, 2)];
  
  ML <- 0;
  ML <- ML + sum(lgamma(As + Bs + 1) - lgamma(As + 1) - lgamma(Bs + 1));
  ML <- ML - sum(lgamma(alpha + beta + As + Bs) - lgamma(alpha + As) - lgamma(beta + Bs));
  ML <- ML +  length(As) * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta));
  
  ML <- ML - lambda * log(alpha + beta);
  
  return(-ML);  
  
}


my_trans <- function(x = 0) {
  return( max(0, -log10(x)));
}