# frechesTest 0.1.0

## New Features

* Initial release of the `frechesTest` package
* Implements ANOVA-style test for jumps in conditional Fréchet means
* Support for multiple metric spaces:
  * Probability distributions (Wasserstein-2 distance)
  * Covariance/correlation matrices (various metrics)
  * Spherical data (geodesic distance)
  * Network Laplacians (for graph-valued data)
* Automatic bandwidth selection via K-fold cross-validation
* Manual bandwidth specification option
* Comprehensive input validation and error handling
* Extensive test suite covering all metric spaces
* OSQP optimization for constrained projection problems

## Functions

* `frechesTest()`: Main function for testing jumps in conditional Fréchet means

## Dependencies

* Core dependencies: `stats`, `pracma`, `Matrix`, `osqp`, `trust`, `fdadensity`, `frechet`
* Suggested packages: `MASS`, `testthat`, `knitr`, `rmarkdown`

## Documentation

* Complete roxygen2 documentation with examples
* README with installation instructions and quick start guide
* Comprehensive help files for all exported functions

## Performance Optimizations

* Pre-computation of quantile functions for density objects
* OSQP model caching for repeated optimizations
* Efficient sparse matrix operations
* Optimized cross-validation implementation