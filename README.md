# estimateW
This is the development repository of the [`R`](https://www.r-project.org/) package `estimateW`.

## Features

The package provides methods to estimate spatial weight matrices in spatial autoregressive type models.

## Installation

Type into your `R` session:

```r
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github(
  repo = "https://github.com/tkrisztin/estimateW")
```

## Demonstration

```r
# Load the package
library(estimateW)

# Estimate a Bayesian model using covid infections data
Y = as.matrix(covid$infections_pc)
X = model.matrix(~infections_pc_lag + stringency_2weekly + precipProbability + temperatureMax + ISO3 + 0,data = covid)
res = sarw(Y = Y,tt = 19,Z = X,niter = 100,nretain = 50)

# Plot the posterior of the spatial weight matrix
plot(res)
```

## References

Tamás Krisztin & Philipp Piribauer (2022) A Bayesian approach for the estimation of weight matrices in spatial autoregressive models, Spatial Economic Analysis, DOI: [`10.1080/17421772.2022.2095426`](https://doi.org/10.1080/17421772.2022.2095426) 
