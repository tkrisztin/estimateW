## Update to CRAN: estimateW 0.1.0

This is an update to the existing CRAN package `estimateW`.

### Changes in this version:

* Added new dataset: `nuts1growth`.
* Added new function: `semw()` for handling spatial error models.
* Corrected labeling of `sigma` to `sigma2` in results.
* Added option to `sim_dgp()` for supplying an exogenous spatial weight matrix.

### R CMD check results

There were no ERRORs or WARNINGs.  
There was one NOTE:

* `checking for future file timestamps ... NOTE`
  - This is due to system time validation in my local environment and is not related to package functionality.  
  - To my knowledge, this is safe to ignore and does not affect CRAN policy compliance.

The package passes `R CMD check --as-cran` on:
- macOS Sonoma (local), R 4.3.2
- Ubuntu 22.04 (Docker), R 4.3.2

I have verified that the package does not write to the user’s home directory or use default write paths. All examples and functions use safe temporary file paths where applicable.
