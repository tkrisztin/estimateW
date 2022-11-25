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
