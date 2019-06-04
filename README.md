# MOG

R package of Bayesian Clustering incorporating multi-layer overlapping groups (MOGC)

## Required Package

Rcpp, RcppEigen, coda, label.switching

### Installing

In R console

```
library(devtools)
install_github("lizhu06/MOGC")
```


## Running the simulations

```
library(MOGC)
seed <- 1
G <- 8
P <- G*30
U1 <- matrix(0, P, G)
temp <- 0
for(g in 1:G){
  U1[temp+seq(1, 30), g] <- 1
  temp <- temp + 30
}
#U1[1, 2] <- 1
#U1[2, 2] <- 1
rho <- 0
strongSignal <- 0.8
weakSignalRatio <- 0.8
strongSignalPerc <- 0.5
weakSignalPerc <- 0.3
g_index <- c(rep(1,2), rep(0, 6))
Sdata = simuSOGC(seed=seed, K=3,
                   numSamplesPerK=c(40,40,40),
                   U1, g_index,  strongSignal, weakSignalRatio,
                   strongSignalPerc, weakSignalPerc,
                   noiseMean=0, sample_sigma=1,
                   sigma_noise=1, percConfounder=0, rho)

Y <- Sdata$data
K <- 3
res_sogc <- SOGC(Y, U1, K=K, center=FALSE,
                 burnInIter=1000, keepIter=3000, maxIter=10000,
                 print_int=10000, debug=FALSE, init_mu=NULL,
                 init_z=NULL, fix_z=FALSE,
                 fix_mu=FALSE, adj_ls=TRUE,
                 seed=123, BernoulliWeighted=TRUE,
                 MH_ind=0, pi_j_prop_n=10)

library(mclust)
(ARI_sogc <- adjustedRandIndex(res_sogc$label_map,Sdata$label))

```

## Authors

* **Li Zhu** - liz86@pitt.edu


