# tests/testthat/test-covariance_space.R
context("frechesTest - Covariance/Correlation Space")

test_that("Covariance space test (Frobenius) runs for NO JUMP", {
  skip_if_not_installed("MASS")
  set.seed(101)
  n_total_cov_nj <- 300 # Increased N
  p_cov <- 3
  X_scalar_cov_nj <- runif(n_total_cov_nj, 0, 1)
  Y_obj_cov_nojump <- vector("list", n_total_cov_nj)
  
  # True Sigma varies smoothly with X
  Sigma_nojump_fun <- function(x_val) {
    val <- 1 + 0.5 * x_val # Smoothly increasing diagonal
    offdiag <- 0.1 + 0.2 * x_val # Smoothly increasing off-diagonal
    mat <- matrix(c(val, offdiag, 0.05, # Example for p_cov=3
                    offdiag, val*1.2, offdiag*0.5,
                    0.05, offdiag*0.5, val*0.8), p_cov, p_cov, byrow=TRUE)
    return(as.matrix(Matrix::nearPD(mat, ensureSymmetry=TRUE)$mat))
  }
  n_samples_for_cov_matrix_nj <- 50
  
  for(i in 1:n_total_cov_nj) {
    true_sigma <- Sigma_nojump_fun(X_scalar_cov_nj[i])
    sample_data <- MASS::mvrnorm(n_samples_for_cov_matrix_nj, mu = rep(0, p_cov), Sigma = true_sigma)
    Y_obj_cov_nojump[[i]] <- stats::cov(sample_data)
  }
  frechet_options_cov_nj <- list(metric = "frobenius")
  
  result_nj <- frechesTest(
    Y_obj = Y_obj_cov_nojump, X_scalar = X_scalar_cov_nj, c_val = 0.5,
    metric_space_type = "covariance",
    h_frechet = "CV", kernel_frechet_char = "epan", # User-specified h for this test
    frechet_optns = frechet_options_cov_nj
  )
  # ... (expectations for no jump: high p-value) ...
  expect_null(result_nj$error)
  expect_true(result_nj$p_value > 0.1, label = "P-value high for cov no jump")
})

test_that("Correlation space test detects a JUMP with refined data", {
  skip_if_not_installed("MASS")
  set.seed(112)
  n1_corr_j <- 75; n2_corr_j <- 75; n_total_corr_j <- n1_corr_j + n2_corr_j
  p_cov <- 3
  X_scalar_corr_j <- c(runif(n1_corr_j, 0, 0.5 - 0.05), runif(n2_corr_j, 0.5 + 0.05, 1))
  Y_obj_corr_jump <- vector("list", n_total_corr_j)
  
  A1 <- matrix(c(1, 0.1, 0.15, 0.1, 1, 0.2, 0.15, 0.2, 1), p_cov, p_cov, byrow = TRUE)
  A2 <- matrix(c(1, 0.7, 0.65, 0.7, 1, 0.8, 0.65, 0.8, 1), p_cov, p_cov, byrow = TRUE)
  n_samples_for_corr_matrix_j <- 150 
  
  for(i in 1:n_total_corr_j) {
    true_underlying_corr <- if (X_scalar_corr_j[i] < 0.5) A1 else A2
    sample_data <- MASS::mvrnorm(n_samples_for_corr_matrix_j, mu = rep(0, p_cov), Sigma = true_underlying_corr)
    Y_obj_corr_jump[[i]] <- stats::cov2cor(stats::cov(sample_data))
  }
  frechet_options_corr_j <- list(metric = "frobenius") 
  
  result_corr_j <- frechesTest(
    Y_obj = Y_obj_corr_jump, X_scalar = X_scalar_corr_j, c_val = 0.5,
    metric_space_type = "correlation",
    h_frechet = "CV", kernel_frechet_char = "epan",
    frechet_optns = frechet_options_corr_j
  )
  # ... (expectations for jump: low p-value) ...
  expect_null(result_corr_j$error)
  expect_true(result_corr_j$p_value < 0.05, label = "P-value low for corr jump")
})

# Log-Cholesky test can remain similar, ensuring Y_obj has variability based on X
test_that("Covariance space Log-Cholesky metric runs with varying objects", {
  skip_if_not_installed("MASS")
  set.seed(131)
  n_total_lc <- 50 # Increased slightly
  p_cov_lc <- 2
  X_scalar_lc <- runif(n_total_lc, 0, 1)
  Y_obj_cov_lc <- vector("list", n_total_lc)
  n_samples_lc <- 40
  for (i in 1:n_total_lc){
    # Sigma's diagonal elements vary with X_scalar_lc
    current_sigma <- diag(c(1 + X_scalar_lc[i]*0.5, 1.5 - X_scalar_lc[i]*0.3))
    # Ensure positive definite for safety, though diag usually is
    current_sigma <- as.matrix(Matrix::nearPD(current_sigma)$mat) 
    sample_data <- MASS::mvrnorm(n_samples_lc, mu = rep(0, p_cov_lc), Sigma = current_sigma)
    Y_obj_cov_lc[[i]] <- stats::cov(sample_data)
  }
  frechet_options_lc <- list(metric = "log_cholesky")
  
  result_lc <- frechesTest(
    Y_obj = Y_obj_cov_lc, X_scalar = X_scalar_lc, c_val = 0.5,
    metric_space_type = "covariance",
    h_frechet = 0.2, kernel_frechet_char = "gauss",
    frechet_optns = frechet_options_lc
  )
  expect_null(result_lc$error)
  expect_true(is.numeric(result_lc$Tn))
})
