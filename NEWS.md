# estimateW 0.2.0

* Added `slx()` function for estimating Spatial Lag of X (SLX) models with exogenous spatial weight matrices.
* Added `logdetSpline()` function implementing spline-based approximations for the log-determinant in rho sampling.
* Updated default rho prior to utilize the Pace and Barry approximation for improved computational efficiency.
* Fixed bug where `print` and `plot` methods failed to correctly handle `SEM` and `SLX` model objects.
* Corrected temporal indices in the `NUTS1` growth dataset.
* Optimized posterior sampling routines for improved convergence in high-dimensional settings.

# v0.1.0 estimateW

- Added new dataset: `nuts1growth`.
- Added new function: `semw()` and `sdmw()`  for handling spatial error models.
- Corrected labeling of `sigma` to `sigma2` in results.
- Added option to `sim_dgp()` for supplying an exogenous spatial weight matrix (defaults to `NULL`).

# v0.0.1 estimateW

- First CRAN release version.
-   Prepare for release with basic functionality
-   `R CMD check --as-cran`: No errors or warnings, one note (New submission)
