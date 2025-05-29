# tests/testthat/test-density_space.R
context("frechesTest - Density Space (Wasserstein)")



test_that("Density space test runs and returns expected structure (no jump)", {
  set.seed(123)
  n_total_ds_nj <- 200 # Moderate N
  X_scalar_ds_nj <- runif(n_total_ds_nj, 0, 1)
  n_obs_per_dist_ds_nj <- 100 # Size of each Y_obj sample
  
  Y_obj_density_nojump <- lapply(1:n_total_ds_nj, function(i) {
    # Continuous mean function of X_scalar
    true_mean <- 0.5 * X_scalar_ds_nj[i] + 0.1 * sin(2 * pi * X_scalar_ds_nj[i])
    true_sd <- 0.3 + 0.1 * X_scalar_ds_nj[i] # SD can also vary
    rnorm(n_obs_per_dist_ds_nj, mean = true_mean, sd = true_sd)
  })
  frechet_options_density_nj <- list(
    qSup = seq(0, 1, length.out = 31),
    den_opts_for_create_density = list(kernelDen = "gauss") 
  )
  
  result_nj <- frechesTest(
    Y_obj = Y_obj_density_nojump, X_scalar = X_scalar_ds_nj, c_val = 0.5,
    metric_space_type = "density", h_frechet = 0.1, kernel_frechet_char = "epan",
    frechet_optns = frechet_options_density_nj, h_fx = 0.1
  )
  expect_null(result_nj$error)
  expect_true(result_nj$p_value > 0.1, label = "P-value high for density no jump")
})

test_that("Density space test detects a clear jump", {
  set.seed(456)
  n1_ds_j <- 100; n2_ds_j <- 100; n_total_ds_j <- n1_ds_j + n2_ds_j
  X_scalar_ds_j <- c(runif(n1_ds_j, 0, 0.48), runif(n2_ds_j, 0.52, 1)) 
  n_obs_per_dist_ds_j <- 40
  
  Y_obj_density_jump <- lapply(1:n_total_ds_j, function(i) {
    true_mean <- if (X_scalar_ds_j[i] < 0.5) 0.0 else 0.8 # Clear jump
    true_sd <- 0.4
    rnorm(n_obs_per_dist_ds_j, mean = true_mean, sd = true_sd)
  })
  frechet_options_density_j <- list(
    qSup = seq(0, 1, length.out = 31),
    den_opts_for_create_density = list(kernelDen = "gauss")
  )
  
  result_j <- frechesTest(
    Y_obj = Y_obj_density_jump, X_scalar = X_scalar_ds_j, c_val = 0.5,
    metric_space_type = "density", h_frechet = 0.11, kernel_frechet_char = "epan",
    frechet_optns = frechet_options_density_j, h_fx = 0.1
  )
  # ... (expectations for jump) ...
  expect_null(result_j$error)
  expect_true(result_j$p_value < 0.05, label = "P-value low for density jump")
})

test_that("Density space handles observation as quantile list", {
  set.seed(789)
  n_total <- 20
  X_scalar <- runif(n_total, 0, 1)
  qSup_orig <- seq(0.05, 0.95, length.out = 21)
  Y_obj_qlist <- lapply(1:n_total, function(i) {
    list(q = stats::qnorm(qSup_orig, mean = X_scalar[i], sd = 1), qSup = qSup_orig, type = "quantile")
  })
  
  frechet_options_density <- list(
    qSup = seq(0, 1, length.out = 31) # Target qSup for test
  )
  
  result <- frechesTest(
    Y_obj = Y_obj_qlist, X_scalar = X_scalar, c_val = 0.5,
    metric_space_type = "density",
    h_frechet = 0.2, kernel_frechet_char = "gauss",
    frechet_optns = frechet_options_density
  )
  expect_null(result$error)
  expect_true(is.numeric(result$p_value))
})


