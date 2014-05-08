library(VGAM);

marginalLikelihood <- function(params, data) {
  
  alpha <- params[1];
  beta <- params[2];
  vec <- data;
  
  vdim <- vec[1];
  As <- vec[1:vdim + 1];
  Bs <- vec[(vdim + 1):(2 * vdim) + 1];
  Ws <- vec[(2 * vdim + 1): (3 * vdim) + 1];
  
  
  ML <- 0;
  ML <- ML + sum( Ws * (lgamma(As + Bs + 1) - lgamma(As + 1) - lgamma(Bs + 1)) );
  ML <- ML - sum( Ws * (lgamma(alpha + beta + As + Bs) - lgamma(alpha + As) - lgamma(beta + Bs)));
  ML <- ML +  sum(Ws) * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta));
  
  ML <- ML - 0.5 * log(alpha + beta);
  
  return(-ML);  
  
}



EMforBetaBinomial <- function(ref, alt) {
  
  param_alpha1 <- 2;
  param_alpha2 <- 10;
  param_pi <- 0.1;
  for (n in 1:1000) {
  
    A <- dbetabinom.ab(alt, ref + alt, param_alpha1, param_alpha2);
    B <- dbetabinom.ab(alt, ref + alt, 30, 30);
    weight <- (1 - param_pi) * A / ((1 - param_pi) * A + param_pi * B);
  
    res <- constrOptim(c(param_alpha1,  param_alpha2), marginalLikelihood, grad=NULL, ui=matrix(c(1, 0, 1, 0, 1, 1), 3, 2), ci=c(0.1, 1, 1), data=c(length(alt), alt, ref, weight));
  
    pre_param_pi <- param_pi;
    param_pi <- 1 - sum(weight) / length(weight);
  
    pre_param_alpha1 <- param_alpha1;
    pre_param_alpha2 <- param_alpha2;
    param_alpha1 <- res$par[1];
    param_alpha2 <- res$par[2];
  
    edelta <- (pre_param_pi - param_pi)^2 + (pre_param_alpha1 - param_alpha1)^2 + (pre_param_alpha2 - param_alpha2)^2;
    # print(c(param_alpha1, param_alpha2, param_pi, edelta));
    if (edelta < 10e-10) {
      break;
    }
  }
  return(c(param_alpha1, param_alpha2, param_pi));
}



generateSimuData <- function(N = 20, D = 50, pi_0 = 0.3, alpha_0 = 1, beta_0 = 10) {
  
  N <- 20;
  D <- 50;
  param_pi_true <- 0.3;
  ref <- rep(0, N);
  alt <- rep(0, N);
  
  temp1 <- rbetabinom.ab(N, D, 50, 50);
  temp2 <- rbetabinom.ab(N, D, alpha_0, beta_0);
  
  Us <- runif(N, 0, 1) ;
  alt[Us < pi_0] <- temp1[Us < pi_0];
  ref[Us < pi_0] <- D - temp1[Us < pi_0];
  alt[Us >= pi_0] <- temp2[Us >= pi_0];
  ref[Us >= pi_0] <- D - temp2[Us >= pi_0];

  return(list(ref, alt));
  
}


evalTimeBetaBinomEst <- function(N = 20, D = 50, param_pi = 0.3, param_alpha = 0.1, param_beta = 10) {

  # N: The number of non-matched references
  # D: Sequencing depth
  # param_pi: The ratio of references having hetero SNPs
  # param_alpha: The first shape parameter in the Beta-Binomial distribution
  # param_beta: The second shape parameter in the Beta-Binomial distribution
  
  TIMEs <- matrix(0, 1000, 3);

for(i in 1:1000) {
  # generating the simulation data
  res <- generateSimuData(N, D, param_pi, param_alpha, param_beta);
  ref <- res[[1]];
  alt <- res[[2]];
  
  # esitimating the parameters using EM algorithm and measuring the computational time
  st <- proc.time();
  params <- EMforBetaBinomial(ref, alt);
  et <- proc.time();
  
  TIMEs[i,] <- as.numeric(et - st)[1:3]
}

hist(TIMEs[,3], breaks=seq(min(TIMEs[,3]) - 0.01, max(TIMEs[,3]) + 0.01, 0.01), col="lightblue", main="", xlab="time");
}

