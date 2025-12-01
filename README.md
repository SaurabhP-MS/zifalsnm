# zifalsnm
**Z**ero **I**nflated **F**actor **A**nalysis **L**ogistic **S**kew **N**ormal **M**ultinomial (*zifalsnm*) model's objective is to achieve effective dimension reduction to manage high dimensionality, accounting for the microbiome data compositional nature, zero inflation and, it utilizes skew-normal priors on the latent factors to explicitly model skewness in log-ratio transformation.

# Installation
To install *zifalsnm* you need to first install *devtool* package. Run the following R code.

```r 
if (!require(devtools)) {
    install.packages("devtools")
    library(devtools)
}
install_github("SaurabhP-MS/zifalsnm", dependencies=TRUE)
library(zifalsnm)
```

# Implementation
Below is an example on how to implement it.

```r
install.packages("sn")
library(sn)

N <- 50 # Number of Samples
P <- 50 # Number of Taxa
K <- 2  # Number of Factors

# Latent Factors
alpha_prior <- 2
alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
f <- matrix(sn::rsn(N * K, alpha_star, 1, alpha_prior), N, K, TRUE)

# Factor Loading
beta <- matrix(runif(P * K, -2.5, 2.5), nrow = P, ncol = K, byrow = TRUE)

kappa <- runif(P)

beta_0 <- runif(P, -2, 2)

z <- matrix(0, N, P, byrow = TRUE)
for (i in 1:N) {
  z[i, ] <- rbinom(P, size = 1, prob = kappa)
}
linear_predictor <- matrix(beta_0, N, P, byrow = TRUE) + tcrossprod(f, beta)
rho <- (1 - z) * exp(linear_predictor) / (rowSums((1 - z) * exp(linear_predictor)) + 1e-8)

cmat <- matrix(0, N, P)

for (i in 1:N) cmat[i, ] <- rmultinom(1, as.integer(runif(1, 80000, 100000)), prob = rho[i, ])

OUTPUT <- ZIFA_LSNM(X = cmat, num_fac = K, NU1P = 1, NU2P = 3, G1P = 200, G2P = 1 / 200, AlphaP = 2)
```




