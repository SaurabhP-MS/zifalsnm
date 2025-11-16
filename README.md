# zifalsnm
**Z**ero **I**nflate **F**actor **A**nalysis **L**ogistic **S**kew **N**ormal **M**ultinomial (*zifalsnm*) model's objective is to achieve effective dimension reduction to manage high dimensionality, accounting for the microbiome data compositional nature, zero inflation and, it utilizes skew-normal priors on the latent factors to explicitly model skewness in log-ratio transformation.

# Installation
To install *zifalsnm* you need to first install *devtool* package. The steps are given below.

```{r setup}
if (!require(devtools)) {
    install.packages("devtools")
    library(devtools)
}
install_github("SaurabhP-MS/zifalsnm", dependencies=TRUE)
library(zifalsnm)
```

# Implementation


