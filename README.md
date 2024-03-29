# estimateW
<!-- badges: start -->
[![R-CMD-check](https://github.com/tkrisztin/estimateW/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/tkrisztin/estimateW/actions/workflows/check-standard.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/estimateW)](https://CRAN.R-project.org/package=estimateW)
[![month](https://cranlogs.r-pkg.org/badges/estimateW)](https://www.r-pkg.org/pkg/estimateW)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/estimateW)](https://www.r-pkg.org/pkg/estimateW)
<!-- badges: end -->

This is the development repository of the [`R`](https://www.r-project.org/) package [`estimateW`](https://cran.r-project.org/package=estimateW).

## Features

The package provides methods to estimate spatial weight matrices in spatial autoregressive type models.

## Install CRAN Version

Type into your `R` session:

```r
install.packages("estimateW")
```

For more information, please visit the [CRAN page](https://cran.r-project.org/package=estimateW) of the package.

## Install Latest Development Version

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
require(dplyr)

tt = length(unique(covid$date))
n = length(unique(covid$ISO3))

# reorder by date and longitude
covid = covid %>% 
  arrange(date, LON) %>%
  mutate(date = as.factor(date))
  
# Benchmark specification from Krisztin and Piribauer (2022) SEA
Y = as.matrix(covid$infections_pc - covid$infections_pc_lag)
X = model.matrix(~infections_pc_lag + stringency_2weekly + 
                   precipProbability + temperatureMax + ISO3 + as.factor(date) + 0,data = covid)

# use a flat prior for W
flat_W_prior = W_priors(n = n,nr_neighbors_prior = rep(1/n,n))

# Estimate a Bayesian model using covid infections data
res = sarw(Y = Y,tt = tt,Z = X,niter = 200,nretain = 50,
           W_prior = flat_W_prior)
           
# Plot the posterior of the spatial weight matrix
dimnames(res$postw)[[2]] = dimnames(res$postw)[[1]] = covid$ISO3[1:n]
plot(res,font=3,cex.axis=0.75,las=2)
```

## References

Tamás Krisztin & Philipp Piribauer (2022) A Bayesian approach for the estimation of weight matrices in spatial autoregressive models, Spatial Economic Analysis, DOI: [`10.1080/17421772.2022.2095426`](https://doi.org/10.1080/17421772.2022.2095426) 
