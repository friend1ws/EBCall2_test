# Patemeter Estimation of Hetero SNP

## Goal

* Estimating the randomness of the ratios of variant short reads at hetero SNP so that we
discriminate somatic mutations from germline porimorphisms efficiently.

## Method

* Here, we assume that the distribution of the variant short reads follows Beta-Binomial distributions. 
* Genotypes are determined by SNP arrays (Affimetrix Human Mapping 250K Nsp Array).
* 13 exome sequences.

## Results

### Exome

#### Normal




parameter
$\alpha = 113.4064, \beta = 122.8334$

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


#### Tumor




parameter
$\alpha = 20.6556, \beta = 22.2761$

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


### Genome

#### Normal




parameter
$\alpha = 61.8145, \beta = 62.2975$

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


#### Tumor




parameter
$\alpha = 18.4874, \beta = 18.6416$

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



### Observation

* The modes in exome data were closer to center than those in genome data.
* The variances in genome data were larger than those in exome data.




parameter
$\alpha = 61.2842, \beta = 61.8634$

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

