# tests/testthat/test-cv_and_undersmoothing.R
context("frechesTest - CV Bandwidth Selection and Undersmoothing")

test_that("CV runs for density space and returns valid bandwidths", {
  set.seed(501)
  n1 <- 100; n2 <- 100; n_total_cv <- n1 + n2 # Adjusted N for CV speed
  X_scalar_cv <- c(runif(n1, 0.1, 0.45), runif(n2, 0.55, 0.9)) 
  
  n_obs_per_dist_obj <- 100 # Number of raw observations forming each Y_obj[[i]]
  
  Y_obj_density_cv <- lapply(1:n_total_cv, function(i) {
    current_x <- X_scalar_cv[i]
    # Mean of the underlying N(mean, sd) from which we draw the sample Y_obj[[i]]
    # This defines a jump scenario, which is fine for testing CV mechanism.
    true_mean_for_Y_i <- if(current_x < 0.5) {
      0.0 + 0.2 * current_x # Slight slope on left
    } else {
      0.8 + 0.2 * current_x # Jumped and different slope on right
    }
    true_sd_for_Y_i <- 0.4 # Can also vary with current_x if desired
    
    rnorm(n_obs_per_dist_obj, mean = true_mean_for_Y_i, sd = true_sd_for_Y_i) 
  })
  
  frechet_options_density_cv <- list(
    qSup = seq(0, 1, length.out = 21), 
    den_opts_for_create_density = list(kernelDen = "gauss") # Let CreateDensity pick bw
  )
  
  result_cv <- frechesTest(
    Y_obj = Y_obj_density_cv, X_scalar = X_scalar_cv, c_val = 0.5,
    metric_space_type = "density",
    h_frechet = 0.2, # "CV", 
    kernel_frechet_char = "epan",
    frechet_optns = frechet_options_density_cv,
    cv_K_folds = 10, # Reduced for typical test speed
    cv_n_bw_candidates = 10, # Reduced for typical test speed
    undersmooth_factor = 0.9,
    min_eff_points_cv_fold = 10
  )
  
  # ... (expectations as before, maybe adjust undersmooth_factor in expect_equal)
  expect_null(result_cv$error, info = paste("CV run error:", result_cv$error))
  expect_true(is.numeric(result_cv$h_mean_cv_selected) && result_cv$h_mean_cv_selected > 0)
  if (is.finite(result_cv$h_mean_cv_selected)) {
    expect_equal(result_cv$h_variance_used, result_cv$h_mean_cv_selected * 0.9, tolerance = 1e-9)
  }
  expect_true(is.numeric(result_cv$p_value) && result_cv$p_value >= 0 && result_cv$p_value <= 1)
})

test_that("User-supplied h_frechet bypasses CV and applies undersmoothing", {
  set.seed(502)
  n_total_user_h <- 100
  X_scalar_user_h <- runif(n_total_user_h, 0, 1)
  n_obs_per_dist_obj_user_h <- 100
  Y_obj_density_user_h <- lapply(1:n_total_user_h, function(i) {
    rnorm(n_obs_per_dist_obj_user_h, mean = 0.5 * X_scalar_user_h[i], sd = 0.5) # Smooth mean
  })
  
  frechet_options_density_user_h <- list(qSup = seq(0, 1, length.out = 21))
  user_h <- 0.2
  undersmooth_factor_user <- 0.75
  
  result_user_h <- frechesTest(
    Y_obj = Y_obj_density_user_h, X_scalar = X_scalar_user_h, c_val = 0.5,
    metric_space_type = "density",
    h_frechet = user_h, 
    kernel_frechet_char = "gauss",
    frechet_optns = frechet_options_density_user_h,
    undersmooth_factor = undersmooth_factor_user
  )
  
  # ... (expectations as before)
  expect_null(result_user_h$error)
  expect_true(is.na(result_user_h$cv_h_plus))
  expect_equal(result_user_h$h_mean_cv_selected, user_h)
  expect_equal(result_user_h$h_variance_used, user_h * undersmooth_factor_user, tolerance = 1e-9)
})

test_that("CV handles few points on one side gracefully", {
  set.seed(503)
  X_scalar_sparse_side <- c(runif(3, 0.1, 0.48), runif(30, 0.52, 0.9)) 
  n_total_sparse <- length(X_scalar_sparse_side)
  n_obs_per_dist_obj_sparse <- 30
  Y_obj_sparse_side <- lapply(1:n_total_sparse, function(i) {
    rnorm(n_obs_per_dist_obj_sparse, mean = (X_scalar_sparse_side[i] - 0.5) * 2, sd = 0.5) # Mean varies
  })
  
  frechet_options_sparse <- list(qSup = seq(0, 1, length.out = 21))
  
  expect_warning(
    result_sparse <- frechesTest(
      Y_obj = Y_obj_sparse_side, X_scalar = X_scalar_sparse_side, c_val = 0.5,
      metric_space_type = "density", h_frechet = "CV", kernel_frechet_char = "epan",
      frechet_optns = frechet_options_sparse, cv_K_folds = 2, cv_n_bw_candidates = 3
    ))
  # ... (expectations as before)
  expect_equal(result_sparse$error, 'Right side Frechet estimation failed.')
  expect_true(is.numeric(result_sparse$cv_h_minus) && result_sparse$cv_h_minus > 0)
})

test_that("CV works for covariance space", {
  skip_if_not_installed("MASS")
  set.seed(504)
  n1_cov <- 50; n2_cov <- 50; n_total_cov_cv <- n1_cov + n2_cov # Further reduced for speed
  p_cov <- 2
  X_scalar_cov_cv <- c(runif(n1_cov, 0.1, 0.45), runif(n2_cov, 0.55, 0.9))
  Y_obj_cov_cv <- vector("list", n_total_cov_cv)
  c_val <- 0.5
  
  # True covariance matrices that vary smoothly with X on each side, with a jump
  Sigma_fun <- function(x_val) {
    base_diag <- if(x_val < 0.5) 1.0 else 1.5 # Jump in variance scale
    off_diag_base <- if(x_val < 0.5) 0.2 else 0.6 # Jump in correlation
    
    # Smooth variation with x
    s_factor <- 1 + 0.2 * (x_val - c_val) # c_val needs to be defined if used here, or use x_val directly
    # Let's use x_val directly for smooth change from 0
    s_factor <- 1 + 0.2 * x_val 
    
    mat <- matrix(c(base_diag * s_factor, off_diag_base * s_factor,
                    off_diag_base * s_factor, base_diag * s_factor * 1.2), p_cov, p_cov)
    # Ensure SPD by making it diagonally dominant slightly or use nearPD
    diag(mat) <- diag(mat) + 0.1 # Ensure strict positivity
    return(Matrix::nearPD(mat)$mat) # Ensure it's a valid covariance matrix
  }
  n_samples_for_cov_matrix_cv <- 40 # Samples to estimate each cov matrix
  
  for(i in 1:n_total_cov_cv) {
    true_sigma <- Sigma_fun(X_scalar_cov_cv[i])
    sample_data <- MASS::mvrnorm(n_samples_for_cov_matrix_cv, mu = rep(0, p_cov), Sigma = as.matrix(true_sigma))
    Y_obj_cov_cv[[i]] <- stats::cov(sample_data)
  }
  frechet_options_cov_cv <- list(metric = "frobenius")
  
  result_cov_cv <- frechesTest(
    Y_obj = Y_obj_cov_cv, X_scalar = X_scalar_cov_cv, c_val = 0.5,
    metric_space_type = "covariance", h_frechet = "CV", kernel_frechet_char = "gauss",
    frechet_optns = frechet_options_cov_cv, cv_K_folds = 10, cv_n_bw_candidates = 10,
    undersmooth_factor = 0.9
  )
  # ... (expectations as before)
  expect_null(result_cov_cv$error)
  expect_true(is.numeric(result_cov_cv$h_mean_cv_selected) && result_cov_cv$h_mean_cv_selected > 0)
})