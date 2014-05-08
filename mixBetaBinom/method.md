# On parameter estimation of mixture of Beta-Binomial Distribution

## Goal

For performing fast somatic mutation calling from whole genome data,
The bottle necks will be

1. Pursing the sequence BAM files and get the candidates for somatic mutation roughly.
2. Estimate the posterior probability of being somatic mutation for each candidate.

The first step takes about several hours for whole genome sequence data 
(If we use APIs in samtools, we have to confirm that..).
For the second step, we will need to estimate parameters of a mixture of Beta-Binomial distributions


If the number of candidates remaining after the first step is 20 thousand 
and the estimation time of Beta-Binomial parameters is 1 second,
the total computational time will be around 60 hours...
This is unbearable.
We should aim for at least 0.01 second for one estimation...


## Experiment using R





### Basic code:


```r
evalTimeBetaBinomEst <- function(N = 20, D = 50, param_pi = 0.3, param_alpha = 0.1, param_beta = 10)
  # N: The number of non-matched references
  # D: Sequencing depth
  # param_pi: The ratio of references having hetero SNPs
  # param_alpha: The first shape parameter in the Beta-Binomial distribution
  # param_beta: The second shape parameter in the Beta-Binomial distribution

  TIMEs <- matrix(0, 100, 3);

  for(i in 1:100) {
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
```



### Non-error prone sites on non-SNP sites : param_pi = 0, param_alpha = 0.2, param_beta = 20;

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


### Non-error prone sites on SNP sites : param_pi = 0.4, param_alpha = 0.2, param_beta = 20;

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 



### Error prone sites on non-SNP sites : param_pi = 0.0, param_alpha = 3, param_beta = 20;

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


### Error prone sites on SNP sites : param_pi = 0.4, param_alpha = 3, param_beta = 20;

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



### Artifacts : param_pi = 0, param_alpha = 5, param_beta = 5;

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


## Method

Basically, the approach for estimation will be 

1. obtain a surrogate function which bounded from above by the likelihood function.
2. maximize the surrogate function.


The following papers seem to be very useful.

* "Estimating a Dirichlet distribution", Minka, from Web, 2000.
* "MM Algorithms for Some Discrete Multivariate Distributions", Zhou et al., JCGS, 2010.


## 
