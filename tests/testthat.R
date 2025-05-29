# tests/testthat.R
library(testthat)
library(frechesTest) # Replace YourPackageName with the actual name if you build it

# Source your R files directly if not building a package for quick testing
# Adjust paths as necessary
# source("../../R/FrechetJumpTest_helpers.R") # Assuming helpers are in the R directory
# source("../../R/FrechetJumpTest_main.R")   # Assuming main function is in the R directory
# Also, ensure the `frechet` package's source files are loaded if they aren't a formal dependency.

test_check("frechesTest") # Or your package name