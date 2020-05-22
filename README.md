# Evolutionary State Space Model
This repository is to implement the approach in the paper ["Evolutionary State-Space Model and Its Application to Time-Frequency Analysis of Local Field Potentials"](http://www3.stat.sinica.edu.tw/preprint/SS-2017-0420_Preprint.pdf)

1. Generate AR(2) parameters corresponding to different frequency bands.

```
phi1 = tran_AR2(1.001, (2/1000)*2*pi) #delta
phi2 = tran_AR2(1.001, (6/1000)*2*pi) #theta
phi3 = tran_AR2(1.001, (8/1000)*2*pi) #alpha
phi4 = tran_AR2(1.001, (15/1000)*2*pi) #beta
phi5 = tran_AR2(1.001, (32/1000)*2*pi) #gamma 
```

2. Generate mixtures of AR(2) time series.

```
set.seed(2016)
n_time <- 1000
trials <- 100
n_channel <- 20
n_indpt <- 3
tau <- rep(1, trials)
sigma <- rep(.1, trials)
module <- evol(c(1.001, 1.001, 1.001), 0.00005, trials, n_indpt)
true_M <- matrix(runif(n_channel*n_indpt, min = 0, max = 1), n_channel, n_indpt)
sim <- generate_time(module, sigma, tau, n_channel, n_indpt, trials, n_time, true_M)
```

3. Simulations 
