## Update to CRAN: estimateW 0.2.0

This is an update to the existing CRAN package `estimateW` (previous version 0.1.0).

### Changes in this version:

* New functionality: Added `slx()` for Spatial Lag of X models and `logdetSplie()` for spline-based log-determinant approximations in rho sampling.
* Algorith update: Updated the default rho prior to utilize the Pace and Barry approximation for improved computational efficiency.
* Bug fixes: Resolved an issue where print and plot methods failed to handle SEM and SLX model objects correctly.
* Data correction: Corrected temporal indices in the `NUTS1` growth dataset.

### R CMD check results

There were no ERRORs or WARNINGs.  
There was one NOTE:

* `checking for future file timestamps ... NOTE`
  - This is due to system time validation in my local environment and is not related to package functionality.  
  - To my knowledge, this is safe to ignore and does not affect CRAN policy compliance.

The package passes `R CMD check --as-cran` on:
- macOS Tahoe (local), R 4.4.3
- Ubuntu 22.04 (Docker), R 4.3.2
- Windows Server 2022 (R-devel and R-release)

I have verified that the package does not write to the user’s home directory or use default write paths. All examples and functions use safe temporary file paths where applicable. No new dependencies have been added that are not available on CRAN.
