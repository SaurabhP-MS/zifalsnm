# zifalsnm
**Z**ero **I**nflated **F**actor **A**nalysis **L**ogistic **S**kew **N**ormal **M**ultinomial (*zifalsnm*) model's objective is to achieve effective dimension reduction to manage high dimensionality, accounting for the microbiome data compositional nature, zero inflation and, it utilizes skew-normal priors on the latent factors to explicitly model skewness in log-ratio transformation. The model is fitted using coordinate ascent variational inference (CAVI).

# Installation
To install *zifalsnm* you need to first install the *remotes* (or *devtools*) package. Run the following R code.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("SaurabhP-MS/zifalsnm", dependencies = TRUE)
library(zifalsnm)
```

# Implementation
Below is an example on how to implement it. The last column of the count matrix `X` is used as the reference taxon. Taxa with only zero counts should be removed before fitting the model.

```r
library(sn)
set.seed(1)

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

# Remove taxa with only zero counts before fitting the model
cmat <- cmat[, colSums(cmat) > 0]

OUTPUT <- ZIFA_LSNM(X = cmat, number_of_factors = K, epsilon = 5e-5,
                    NU1Prior = 1, NU2Prior = 3, G1Prior = 200, G2Prior = 1 / 200,
                    AlphaPrior = 2, Max_Iter = 200, verbose = TRUE)
```

## Arguments

| Argument | Description | Default |
|---|---|---|
| `X` | Count matrix (samples in rows, taxa in columns); the last column is the reference taxon | – |
| `number_of_factors` | Number of latent factors | `2` |
| `epsilon` | Convergence tolerance on the relative change in the ELBO | `5e-5` |
| `NU1Prior`, `NU2Prior` | Prior shape parameters of the zero-inflation probability κ<sub>j</sub> | `1`, `3` |
| `G1Prior`, `G2Prior` | Prior shape and rate parameters of δ<sub>jt</sub> | `200`, `1/200` |
| `AlphaPrior` | Prior skewness parameter of the latent factors | `2` |
| `Max_Iter` | Maximum number of CAVI iterations | `200` |
| `verbose` | Print progress messages at each iteration | `TRUE` |

The output is a list with the variational parameters (`R`, `Lambda`, `Xi`, `Omega`, `Alpha`, `A0`, `C0`, `G1`, `G2`, `Tau_1`, `Tau_2`, `Pi_ij`), the `Estimated_Compositions`, the `ELBO_Trace` and the trace of every parameter across iterations.

# Reproducibility
All code needed to reproduce the results in the manuscript is in the [`Reproducible_R_Codes`](Reproducible_R_Codes) folder. Set your R working directory to that folder before running the scripts.

| File | Reproduces |
|---|---|
| `Reproducible_Code_For_Data_Analysis.R` | Pre-processing of the Jacobs et al. (2016) IBD data, the skewness plot, and the real data analysis (factor loadings and latent factor plots) |
| `Jacobs_ibd_2016` | The publicly available IBD data used in the real data analysis |
| `Reproducible_Code_For_Violin_Plots_k=2.R` | Violin plots of the RMSE values for the simulations with K = 2 |
| `Reproducible_Code_For_Violin_Plots_k=5.R` | Violin plots of the RMSE values for the simulations with K = 5 |
| `RMSE_Values_k=2/`, `RMSE_Values_k=5/` | Saved RMSE values from the simulation study |
