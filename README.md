# zifalsnm
**Z**ero **I**nflated **F**actor **A**nalysis **L**ogistic **S**kew **N**ormal **M**ultinomial (*zifalsnm*) model's objective is to achieve effective dimension reduction to manage high dimensionality, accounting for the microbiome data compositional nature, zero inflation and, it utilizes skew-normal priors on the latent factors to explicitly model skewness in log-ratio transformation. The model is fitted using coordinate ascent variational inference (CAVI).

# Installation
To install *zifalsnm* you need to first install the *remotes* (or *devtools*) package. Run the following R code.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
options(timeout = max(1200, getOption("timeout")))
remotes::install_github("SaurabhP-MS/zifalsnm", dependencies = TRUE)
library(zifalsnm)
```

# Implementation
Below is an example on how to implement it. The last column of the count matrix `X` is used as the reference taxon.

```r

OUTPUT <- ZIFA_LSNM(X = Count_Matrix, number_of_factors = K, epsilon = 5e-5,
                    NU1Prior = 1, NU2Prior = 3, G1Prior = 1.5, G2Prior = 1.5,
                    AlphaPrior = 2, Max_Iter = 200, verbose = TRUE)
```

## Arguments

| Argument | Description | Default |
|---|---|---|
| `X` | Count matrix (samples in rows, taxa in columns); the last column is the reference taxon | – |
| `number_of_factors` | Number of latent factors | `2` |
| `epsilon` | Convergence tolerance on the relative change in the ELBO | `5e-5` |
| `NU1Prior`, `NU2Prior` | Prior shape parameters of the zero-inflation probability κ<sub>j</sub> | `1`, `3` |
| `G1Prior`, `G2Prior` | Prior shape and rate parameters of δ<sub>jt</sub> | `1.5`, `1.5` |
| `AlphaPrior` | Prior skewness parameter of the latent factors | `2` |
| `Max_Iter` | Maximum number of CAVI iterations | `200` |
| `verbose` | Print progress messages at each iteration | `TRUE` |

The output is a list with the variational parameters (`R`, `Lambda`, `Xi`, `Omega`, `Alpha`, `A0`, `C0`, `G1`, `G2`, `Tau_1`, `Tau_2`, `Pi_ij`), the `Estimated_Compositions`, the `ELBO_Trace` and the trace of every parameter across iterations.
