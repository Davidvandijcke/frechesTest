# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R package called `frechesTest` that implements an ANOVA-style test for detecting discontinuities (jumps) in conditional Fréchet means of random objects in metric spaces. The package supports multiple metric spaces including probability distributions (Wasserstein distance), covariance/correlation matrices, spherical data, and network Laplacians.

## Common Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (creates frechesTest_0.1.0.tar.gz)
R CMD check frechesTest_0.1.0.tar.gz

# Install locally
R CMD INSTALL .
```

### Using devtools (preferred for development)
```r
# Load all functions for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Run a specific test file
devtools::test(filter = "covariance")  # runs test-covariance_space.R

# Check the package
devtools::check()

# Build documentation
devtools::document()

# Install the package
devtools::install()
```

### Running Tests
```r
# Run all tests
testthat::test_local()

# Run specific test file
testthat::test_file("tests/testthat/test-general_functionality.R")
```

## Architecture and Key Components

### Core Function Structure
The package revolves around a single main function `frechesTest()` which:
1. Accepts metric-space valued data and a scalar running variable
2. Tests for jumps at a specified cutoff using local polynomial regression
3. Implements bandwidth selection via cross-validation
4. Returns test statistics and p-values

### Key Implementation Details
- **Main function**: `R/frechesTest.R` contains the primary `frechesTest()` function that orchestrates the test
- **Helper functions**: `R/frechesTest_helpers.R` contains internal functions for:
  - One-sided local polynomial weight calculation
  - Fréchet mean and variance estimation
  - Bandwidth selection via K-fold cross-validation
  
### Metric Space Handling
The package uses the `frechet` package for metric space computations. Different metric spaces require different preprocessing:
- **Density space**: Converts raw data to quantile functions
- **Covariance space**: Ensures positive definiteness
- **Sphere space**: Normalizes vectors to unit length
- **Network space**: Handles graph Laplacian matrices

### Testing Structure
Tests are organized by metric space type:
- `test-general_functionality.R`: Core functionality and edge cases
- `test-density_space.R`: Wasserstein distance tests
- `test-covariance_space.R`: Matrix metric tests
- `test-sphere_space.R`: Spherical data tests
- `test-cv_and_undersmoothing.R`: Bandwidth selection tests

Each test file contains multiple test cases verifying both statistical properties (power, size) and implementation correctness.