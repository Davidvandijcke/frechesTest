# tests/testthat/test-general_functionality.R
context("frechesTest - General Functionality and Edge Cases")

test_that("Handles insufficient data on one side", {
  set.seed(301)
  n_total <- 20
  X_scalar_left_sparse <- c(runif(1, 0, 0.49), runif(n_total-1, 0.51, 1)) # Only 1 point on left
  Y_obj_density_edge <- lapply(1:n_total, function(i) rnorm(10, mean=X_scalar_left_sparse[i]))
  
  # Valid qSup for density type
  frechet_options_density_valid_qSup <- list(
    qSup = seq(0, 1, length.out = 21) 
  )
  
  expect_warning(
    result_left_fail <- frechesTest(
      Y_obj = Y_obj_density_edge, X_scalar = X_scalar_left_sparse, c_val = 0.5,
      metric_space_type = "density",
      h_frechet = 0.1, kernel_frechet_char = "epan",
      frechet_optns = frechet_options_density_valid_qSup
    ),
    "Too few points on side - for local polynomial regression"
  )
  expect_true(is.na(result_left_fail$Tn) || !is.null(result_left_fail$error))
  
  X_scalar_right_sparse <- c(runif(n_total-1, 0, 0.49), runif(1, 0.51, 1)) # Only 1 point on right
  # Y_obj needs to match X_scalar_right_sparse length and values for this specific test
  Y_obj_density_edge_rev <- lapply(1:n_total, function(i) rnorm(10, mean=X_scalar_right_sparse[i]))
  
  expect_warning(
    result_right_fail <- frechesTest(
      Y_obj = Y_obj_density_edge_rev, X_scalar = X_scalar_right_sparse, c_val = 0.5,
      metric_space_type = "density",
      h_frechet = 0.1, kernel_frechet_char = "epan",
      frechet_optns = frechet_options_density_valid_qSup
    ),
    "Too few points on side \\+ for local polynomial regression"
  )
  expect_true(is.na(result_right_fail$Tn) || !is.null(result_right_fail$error))
})

test_that("Kernel constant S_K calculation for different kernels", {
  set.seed(311) 
  X_scalar_dummy <- runif(30,0,1) # Increased n for stability
  Y_obj_dummy <- lapply(X_scalar_dummy, function(x) rnorm(20,x)) # Increased obs per Y_obj
  # Valid qSup for density type
  frechet_options_dummy_valid_qSup <- list(qSup = seq(0, 1, length.out=21)) 
  
  # Ensure enough points on both sides for these tests to not fail due to sparse data
  # For S_K test, we just need the function to run through, actual values depend on data.
  # We can use a symmetric X around c_val to ensure both sides are populated.
  X_scalar_sk_test <- c(runif(15, 0.0, 0.45), runif(15, 0.55, 1.0))
  Y_obj_sk_test    <- lapply(X_scalar_sk_test, function(x) rnorm(20, x))
  
  
  res_gauss <- frechesTest(Y_obj_sk_test, X_scalar_sk_test, 0.5, "density", 0.15, "gauss", frechet_options_dummy_valid_qSup)
  expect_null(res_gauss$error, info = "Gaussian kernel S_K test should not error")
  expect_false(is.na(res_gauss$S_K), info = "S_K for Gaussian kernel should not be NA")
  expect_true(is.finite(res_gauss$S_K), info = "S_K for Gaussian kernel should be finite")
  expect_true(res_gauss$S_K > 0, info = "S_K for Gaussian kernel should be positive")
  
  res_epan <- frechesTest(Y_obj_sk_test, X_scalar_sk_test, 0.5, "density", 0.15, "epan", frechet_options_dummy_valid_qSup)
  expect_null(res_epan$error, info = "Epanechnikov kernel S_K test should not error")
  expect_false(is.na(res_epan$S_K), info = "S_K for Epanechnikov kernel should not be NA")
  expect_true(is.finite(res_epan$S_K), info = "S_K for Epanechnikov kernel should be finite")
  expect_true(res_epan$S_K > 0, info = "S_K for Epanechnikov kernel should be positive")
  
  res_rect <- frechesTest(Y_obj_sk_test, X_scalar_sk_test, 0.5, "density", 0.15, "rect", frechet_options_dummy_valid_qSup)
  expect_null(res_rect$error, info = "Rectangular kernel S_K test should not error")
  expect_false(is.na(res_rect$S_K), info = "S_K for Rectangular kernel should not be NA")
  expect_true(is.finite(res_rect$S_K), info = "S_K for Rectangular kernel should be finite")
  expect_true(res_rect$S_K > 0, info = "S_K for Rectangular kernel should be positive")
})

test_that("f_X(c) estimation works", {
  set.seed(321)
  X_scalar_fx <- rnorm(100, mean = 0.5, sd = 0.1) 
  Y_obj_dummy_fx <- lapply(X_scalar_fx, function(x) rnorm(10,x))
  # Valid qSup for density type
  frechet_options_dummy_fx_qSup <- list(qSup = seq(0, 1, length.out=11))
  
  res_fx_default_h <- frechesTest(Y_obj_dummy_fx, X_scalar_fx, 0.5, "density", 0.1, "gauss", frechet_options_dummy_fx_qSup, h_fx = NULL)
  expect_null(res_fx_default_h$error)
  expect_false(is.na(res_fx_default_h$f_X_hat_c))
  expect_true(res_fx_default_h$f_X_hat_c > 0)
  expect_true(res_fx_default_h$f_X_hat_c > 1 && res_fx_default_h$f_X_hat_c < 10, 
              info = paste("f_X_hat_c was", res_fx_default_h$f_X_hat_c, "expected ~3.98"))
  
  res_fx_user_h <- frechesTest(Y_obj_dummy_fx, X_scalar_fx, 0.5, "density", 0.1, "gauss", frechet_options_dummy_fx_qSup, h_fx = 0.05)
  expect_null(res_fx_user_h$error)
  expect_false(is.na(res_fx_user_h$f_X_hat_c))
  expect_true(res_fx_user_h$f_X_hat_c > 0)
  
  X_scalar_edge <- runif(100, 0, 1)
  Y_obj_dummy_edge_fx <- lapply(X_scalar_edge, function(x) rnorm(10,x))

})

test_that("Input validation for qSup in density space works", {
  set.seed(333)
  X_scalar_qSup_test <- runif(20, 0, 1)
  Y_obj_qSup_test <- lapply(X_scalar_qSup_test, function(x) rnorm(10, x))
  
  # qSup not spanning [0,1]
  frechet_optns_bad_qSup_range <- list(qSup = seq(0.1, 0.9, length.out = 11))
  expect_error(
    frechesTest(Y_obj_qSup_test, X_scalar_qSup_test, 0.5, "density", 0.1, "gauss", frechet_optns_bad_qSup_range),
    "frechet_optns\\$qSup must be a numeric vector spanning \\[0, 1]"
  )
  
  # qSup not provided
  frechet_optns_no_qSup <- list()
  expect_error(
    frechesTest(Y_obj_qSup_test, X_scalar_qSup_test, 0.5, "density", 0.1, "gauss", frechet_optns_no_qSup),
    "frechet_optns\\$qSup must be provided"
  )
  

})



test_that("Spherical Geodesic Distance (SpheGeoDist) is correct", {
  # Requires l2norm to be available from frechet package
  l2norm_internal <- frechet:::l2norm 
  
  # Test cases for SpheGeoDist
  p1 <- c(1, 0, 0)
  p2 <- c(0, 1, 0) # Orthogonal, distance should be pi/2
  p3 <- c(-1, 0, 0) # Antipodal, distance should be pi
  p4 <- c(1,0,0) # Identical, distance should be 0
  p5 <- c(1/sqrt(2), 1/sqrt(2), 0) # 45 degrees from p1 in xy plane
  
  expect_equal(frechet::SpheGeoDist(p1, p2), pi/2)
  expect_equal(frechet::SpheGeoDist(p1, p3), pi)
  expect_equal(frechet::SpheGeoDist(p1, p4), 0)
  expect_equal(frechet::SpheGeoDist(p1, p5), pi/4)
  
  # Test normalization robustness
  p6_unnorm <- c(2,0,0)
  p7_unnorm <- c(0,3,0)
  # SpheGeoDist internally normalizes, so result should be pi/2
  expect_equal(frechet::SpheGeoDist(p6_unnorm/l2norm_internal(p6_unnorm), p7_unnorm/l2norm_internal(p7_unnorm)), pi/2)
  # The frechet::SpheGeoDist from the package already has this check and normalization
  # So, a direct call with unnormalized vectors if the package handles it would be:
  # expect_equal(frechet::SpheGeoDist(c(2,0,0), c(0,3,0)), pi/2) # If SpheGeoDist normalizes inputs
  
  # Test for clamping of dot product if numerical precision makes it > 1 or < -1
  # Create two points very close to each other
  slight_pert <- 1e-9
  p_close1 <- c(1,0,0)
  p_close2 <- c(sqrt(1-slight_pert^2), slight_pert, 0)
  p_close2_norm <- p_close2 / l2norm_internal(p_close2)
  dot_prod_close <- sum(p_close1 * p_close2_norm)
  # If dot_prod_close is > 1 due to precision, acos would fail without clamping
  expect_no_error(frechet::SpheGeoDist(p_close1, p_close2_norm)) 
  expect_true(frechet::SpheGeoDist(p_close1, p_close2_norm) >= 0)
  
  # Create two points very nearly antipodal
  p_antip1 <- c(1,0,0)
  p_antip2 <- c(-sqrt(1-slight_pert^2), slight_pert, 0)
  p_antip2_norm <- p_antip2 / l2norm_internal(p_antip2)
  dot_prod_antip <- sum(p_antip1 * p_antip2_norm)
  # If dot_prod_antip is < -1 due to precision, acos would fail without clamping
  expect_no_error(frechet::SpheGeoDist(p_antip1, p_antip2_norm))
  expect_true(frechet::SpheGeoDist(p_antip1, p_antip2_norm) <= pi)
})




test_that("One-sided local polynomial weights sum to 1 over active side", {
  set.seed(123)
  n_test <- 50
  x_obs_test <- runif(n_test, 0, 1)
  c_val_test <- 0.5
  h_val_test <- 0.15
  
  # Assuming K_frechet_fun is correctly loaded and is the actual kernel function
  # and not the character string.
  # K_frechet_fun_epan <- get("kerFctn", envir = asNamespace("frechet"))("epan")
  # K_frechet_fun_gauss <- get("kerFctn", envir = asNamespace("frechet"))("gauss")
  
  # Need to define K_frechet_fun if this test is standalone
  # Sourcing the helpers should make .calculate_one_sided_locpoly_weights available
  # And it internally calls K_frechet_fun which itself calls frechet::kerFctn
  
  # Test for side = "+"
  weights_plus <- .calculate_one_sided_locpoly_weights(
    x_obs = x_obs_test, 
    c_val = c_val_test, 
    h_val = h_val_test, 
    K_fun = frechet:::kerFctn("epan"), # Pass the actual function
    side = "+"
  )
  active_plus_indices <- which(x_obs_test >= c_val_test & weights_plus != 0) # Points contributing to plus-side weights
  # It's possible some points >= c_val still get zero weight if outside kernel support effectively
  # The weights should sum to 1 over the points that actually GOT a non-zero weight from the calculation.
  # However, the standard definition is that they sum to 1 over the points *used in the S_j summations*.
  
  # The weights returned by .calculate_one_sided_locpoly_weights are non-zero only for active_indices_side
  # And they are the w_i^* such that sum(w_i^*) = 1
  if (length(active_plus_indices) >= 2) { # Ensure the weights calculation didn't trivially return all zeros
    expect_equal(sum(weights_plus), 1, tolerance = 1e-9, 
                 label = "Sum of '+' side weights should be 1")
  } else {
    expect_equal(sum(weights_plus), 0, tolerance = 1e-9,
                 label = "Sum of '+' side weights should be 0 if too few points")
  }
  
  # Test for side = "-"
  weights_minus <- .calculate_one_sided_locpoly_weights(
    x_obs = x_obs_test, 
    c_val = c_val_test, 
    h_val = h_val_test, 
    K_fun = frechet:::kerFctn("epan"), 
    side = "-"
  )
  active_minus_indices <- which(x_obs_test < c_val_test & weights_minus != 0)
  if (length(active_minus_indices) >= 2) {
    expect_equal(sum(weights_minus), 1, tolerance = 1e-9,
                 label = "Sum of '-' side weights should be 1")
  } else {
    expect_equal(sum(weights_minus), 0, tolerance = 1e-9,
                 label = "Sum of '-' side weights should be 0 if too few points")
  }
  
  # Test with Gaussian kernel
  weights_plus_gauss <- .calculate_one_sided_locpoly_weights(
    x_obs = x_obs_test, 
    c_val = c_val_test, 
    h_val = h_val_test, 
    K_fun = frechet:::kerFctn("gauss"), 
    side = "+"
  )
  active_plus_indices_gauss <- which(x_obs_test >= c_val_test & weights_plus_gauss != 0)
  if (length(active_plus_indices_gauss) >= 2) {
    expect_equal(sum(weights_plus_gauss), 1, tolerance = 1e-9,
                 label = "Sum of '+' side weights (Gauss) should be 1")
  } else {
    expect_equal(sum(weights_plus_gauss), 0, tolerance = 1e-9,
                 label = "Sum of '+' side weights (Gauss) should be 0 if too few points")
  }
})