## Resubmission
This is a resubmission. In this version I have:

* Converted the DESCRIPTION title to title case.

* Added additional description and a doi reference to the description field in DESCRIPTION 

* Removed changes to the user's par options from plotting

> Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies.
Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir().

* I am aware of this point, and to the best of my knowledge I do not do this: the package does not write, as is CRAN policy. I have checked again, without finding the mentioned problems. I am unsure which line of code gave the wrong impression or what I am overlooking. Please let me know and I will happily fix it.

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* This is a new release.
