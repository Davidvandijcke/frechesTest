# Ensure necessary packages are available if not already loaded by frechet
# library(osqp) # For Wasserstein space mean projection (used internally by frechet)
# library(trust) # For Spherical space mean (used internally by frechet)
# library(Matrix) # For sparse matrices, nearPD (used internally by frechet)
# library(pracma) # For trapz (used internally by frechet)
# library(stats)  # For density, approx, pchisq, integrate, bw.nrd0

# Helper for NULL-coalescing operator (if R version < 4.0.0)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Helper function to get quantiles from various observation types for "density" space
# obs: a single observation (e.g. raw data vector, or list {q, qSup, type="quantile"})
# qSup_target: the common qSup grid to project/calculate quantiles on
# den_opts: options for CreateDensity (like kernelDen, userBwMu, etc.)
.get_quantiles_for_obs_jump_test <- function(obs, qSup_target, den_opts) {
  if (is.list(obs) && !is.null(obs$type) && obs$type == "quantile" && !is.null(obs$q) && !is.null(obs$qSup)) {
    if (length(obs$qSup) == length(qSup_target) && all(obs$qSup == qSup_target)) {
      return(obs$q)
    } else {
      # Interpolate to common qSup_target
      return(stats::approx(x = obs$qSup, y = obs$q, xout = qSup_target, rule = 2)$y)
    }
  } else if (is.numeric(obs)) { # Assumed raw data vector
    current_den_opts <- den_opts %||% list()
    if (!is.null(current_den_opts$kernelDen)) { # Align opt names if needed
      names(current_den_opts)[names(current_den_opts) == "kernelDen"] <- "kernel"
    }
    if (!is.null(current_den_opts$bwDen)) {
      names(current_den_opts)[names(current_den_opts) == "bwDen"] <- "userBwMu"
    }
    
    density_obj <- frechet::CreateDensity(y = obs, optns = current_den_opts)
    quantiles <- fdadensity::dens2quantile(dens = density_obj$y, dSup = density_obj$x, qSup = qSup_target)
    return(quantiles)
  } else {
    stop("Unsupported observation type for .get_quantiles_for_obs_jump_test.")
  }
}


# Helper function to calculate one-sided local polynomial weights
.calculate_one_sided_locpoly_weights <- function(x_obs, c_val, h_val, K_fun, side, n_total) {
  Kh_X_minus_c <- K_fun((x_obs - c_val) / h_val) / h_val
  
  mu_hat_j <- numeric(3) # for j = 0, 1, 2
  indicator <- if (side == "+") {
    x_obs >= c_val
  } else {
    x_obs < c_val
  }
  
  # Calculate mu_hat_j only over the indicated side
  # The n_total in denominator is from paper's definition of hat_L having 1/n
  # If mu_hat_j are defined as sum (not average), then s_in doesn't have n_total in sigma_0_sq_hat
  # Let's follow Fan & Gijbels / Petersen & Mueller (2019) more closely for local poly weights definition:
  # mu_hat_j = sum_{k on side} K_h(X_k-c)(X_k-c)^j
  # Then s_in is defined using these summed mu_hat_j.
  # The 1/n factor is applied when computing the weighted sum for L_hat.
  
  active_indices_side <- which(indicator)
  if (length(active_indices_side) < 2) { # Need at least 2 points for local linear
    warning(paste("Too few points on side", side, "for local polynomial regression. Weights will be zero."))
    return(rep(0, length(x_obs)))
  }
  
  x_obs_side <- x_obs[active_indices_side]
  Kh_X_minus_c_side <- K_fun((x_obs_side - c_val) / h_val) / h_val # K_h values for points on the side
  
  for (j_idx in 0:2) {
    mu_hat_j[j_idx + 1] <- sum(Kh_X_minus_c_side * (x_obs_side - c_val)^j_idx)
  }
  
  sigma_0_sq_hat <- mu_hat_j[1] * mu_hat_j[3] - mu_hat_j[2]^2
  
  s_in <- rep(0, length(x_obs))
  if (abs(sigma_0_sq_hat) < 1e-10) {
    warning(paste("sigma_0_sq_hat is close to zero for side", side, ". Weights might be unstable or zero."))
    # Potentially return uniform weights over active_indices_side if desired, or just zero.
    # For now, this will lead to zero weights if sigma_0_sq_hat is effectively zero.
  } else {
    s_in[active_indices_side] <- (Kh_X_minus_c_side / sigma_0_sq_hat) * (mu_hat_j[3] - mu_hat_j[2] * (x_obs_side - c_val))
  }
  return(s_in)
}

# Helper function to estimate weighted Fréchet mean (generalized for different metric spaces)
.estimate_weighted_frechet_mean <- function(Y_obj_list, active_idx, s_in_weights_active, 
                                            metric_space_type, frechet_optns) {
  if (length(active_idx) == 0) {
    warning("No active samples for weighted Frechet mean estimation.")
    return(NULL)
  }
  Y_active <- Y_obj_list[active_idx]
  
  l_hat_mean <- NULL
  sum_s_in_active <- sum(s_in_weights_active)
  
  if (abs(sum_s_in_active) < 1e-9 && metric_space_type != "sphere") { # For sphere, trust handles weights
    warning("Sum of active weights is close to zero. Fréchet mean may be unstable. Using unweighted mean of active samples.")
    # Fallback to unweighted mean of active objects
    if (metric_space_type == "density") {
      qSup <- frechet_optns$qSup %||% stop("qSup missing for density fallback mean")
      den_opts_cfd <- frechet_optns$den_opts_for_create_density %||% list()
      qin_active_fallback <- t(sapply(Y_active, .get_quantiles_for_obs_jump_test, qSup_target = qSup, den_opts = den_opts_cfd))
      l_hat_q_unprojected <- colMeans(qin_active_fallback)
    } else if (metric_space_type %in% c("covariance", "correlation")) {
      l_hat_mat_unproj <- Reduce("+", Y_active) / length(Y_active)
    } else { # sphere will be handled by trust with potentially zero sum of weights if not careful
      stop("Fallback for zero sum of weights not fully implemented for this space.")
    }
  } else { # Normal case or sphere (trust handles weights)
    if (metric_space_type == "density") {
      qSup <- frechet_optns$qSup %||% stop("qSup missing for density mean")
      den_opts_cfd <- frechet_optns$den_opts_for_create_density %||% list()
      qin_active <- t(sapply(Y_active, .get_quantiles_for_obs_jump_test, qSup_target = qSup, den_opts = den_opts_cfd))
      
      l_hat_q_unprojected <- colSums(qin_active * s_in_weights_active) / sum_s_in_active
      
      m_dim <- ncol(qin_active); Pmat_qp <- as(diag(m_dim), "sparseMatrix")
      A_mono_vals <- rep(c(-1, 1), m_dim - 1); A_mono_i <- rep(1:(m_dim - 1), each = 2)
      A_mono_j <- unlist(lapply(1:(m_dim - 1), function(k) c(k, k + 1)))
      Amat_qp_mono <- Matrix::sparseMatrix(i = A_mono_i, j = A_mono_j, x = A_mono_vals, dims = c(m_dim - 1, m_dim))
      l_bounds_mono <- rep(0, m_dim - 1); u_bounds_mono <- rep(Inf, m_dim - 1)
      
      A_bounds_user <- NULL; l_bounds_user <- NULL; u_bounds_user <- NULL
      if (!is.null(frechet_optns$lower)) {
        A_bounds_user <- rbind(A_bounds_user, Matrix::sparseMatrix(i=1,j=1,x=1, dims=c(1,m_dim)))
        l_bounds_user <- c(l_bounds_user, frechet_optns$lower); u_bounds_user <- c(u_bounds_user, Inf)
      }
      if (!is.null(frechet_optns$upper)) {
        A_bounds_user <- rbind(A_bounds_user, Matrix::sparseMatrix(i=1,j=m_dim,x=1, dims=c(1,m_dim)))
        l_bounds_user <- c(l_bounds_user, -Inf); u_bounds_user <- c(u_bounds_user, frechet_optns$upper)
      }
      Final_Amat_qp <- rbind(Amat_qp_mono, A_bounds_user)
      Final_l_bounds_qp <- c(l_bounds_mono, l_bounds_user); Final_u_bounds_qp <- c(u_bounds_mono, u_bounds_user)
      
      osqp_sol <- osqp::solve_osqp(P=Pmat_qp, q = -l_hat_q_unprojected, A=Final_Amat_qp, l=Final_l_bounds_qp, u=Final_u_bounds_qp, pars = osqp::osqpSettings(verbose = FALSE))
      l_hat_q_projected <- osqp_sol$x
      l_hat_mean <- list(q = l_hat_q_projected, qSup = qSup, type="quantile")
      
    } else if (metric_space_type %in% c("covariance", "correlation")) {
      M_active_list <- Y_active # List of matrices
      metric <- frechet_optns$metric %||% "frobenius"; alpha <- frechet_optns$alpha %||% 1
      
      if (metric %in% c("frobenius", "power")) {
        if (alpha == 0) { # Log-Euclidean like
          log_M_active <- lapply(M_active_list, function(mat) {eig <- eigen(mat, symmetric = TRUE); P <- eig$vectors; Lambda_log <- diag(log(pmax(1e-30, eig$values)), nrow=length(eig$values)); P %*% Lambda_log %*% t(P) })
          mean_log_M <- Reduce("+", mapply("*", log_M_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
          eig_mean_log <- eigen(mean_log_M, symmetric = TRUE); l_hat_mat_unproj <- eig_mean_log$vectors %*% diag(exp(eig_mean_log$values), nrow=length(eig_mean_log$values)) %*% t(eig_mean_log$vectors)
        } else { # Power metric (alpha > 0 includes Frobenius for alpha=1)
          M_alpha_active <- lapply(M_active_list, function(mat) {eig <- eigen(mat, symmetric = TRUE); P <- eig$vectors; Lambda_alpha <- diag(pmax(0, eig$values)^alpha, nrow=length(eig$values)); P %*% Lambda_alpha %*% t(P) })
          mean_M_alpha <- Reduce("+", mapply("*", M_alpha_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
          eig_mean_alpha <- eigen(mean_M_alpha, symmetric = TRUE); l_hat_mat_unproj <- eig_mean_alpha$vectors %*% diag(pmax(0, eig_mean_alpha$values)^(1/alpha), nrow=length(eig_mean_alpha$values)) %*% t(eig_mean_alpha$vectors)
        }
      } else if (metric == "log_cholesky") {
        chol_parts_active <- lapply(M_active_list, function(mat) {LL <- chol(mat); L_part <- LL - diag(diag(LL)); D_diag_part <- diag(LL); list(L = L_part, D_log_diag = log(pmax(1e-30,D_diag_part)))})
        mean_L <- Reduce("+", mapply("*", lapply(chol_parts_active, `[[`, "L"), s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        mean_D_log_diag <- Reduce("+", mapply("*", lapply(chol_parts_active, `[[`, "D_log_diag"), s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        SS <- mean_L + diag(exp(mean_D_log_diag))
        l_hat_mat_unproj <- t(SS) %*% SS
      } else if (metric == "cholesky") {
        L_chol_active <- lapply(M_active_list, chol)
        mean_L_chol <- Reduce("+", mapply("*", L_chol_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        l_hat_mat_unproj <- t(mean_L_chol) %*% mean_L_chol
      } else { stop(paste("Unsupported metric:", metric)) }
    } else if (metric_space_type == "sphere") {
      yin_active_mat <- do.call(rbind, Y_active)
      y0_unnorm <- colSums(yin_active_mat * s_in_weights_active) 
      # If sum_s_in_weights_active is zero, y0_unnorm could be zero or non-zero.
      # If y0_unnorm is zero vector, l2norm is 0.
      norm_y0_unnorm <- frechet:::l2norm(y0_unnorm)
      if (any(is.na(y0_unnorm)) || norm_y0_unnorm < 1e-9) { 
        y0 <- yin_active_mat[1,] # Fallback
      } else {
        y0 <- y0_unnorm / norm_y0_unnorm
      }
      
      objFctn_sphere = function(y_sphere){
        if (!isTRUE(all.equal(frechet:::l2norm(y_sphere),1))) {return(list(value = Inf, gradient=rep(Inf, length(y_sphere)), hessian=diag(Inf, length(y_sphere))))}
        dists_sq_vals <- sapply(1:nrow(yin_active_mat), function(i) frechet::SpheGeoDist(yin_active_mat[i,], y_sphere)^2)
        f_val <- sum(s_in_weights_active * dists_sq_vals)
        
        grads_i_vals <- t(sapply(1:nrow(yin_active_mat), function(i) {
          dist_val <- frechet::SpheGeoDist(yin_active_mat[i,], y_sphere)
          if (abs(dist_val) < 1e-8 || abs(dist_val - pi) < 1e-8) { # Avoid issues at 0 or pi
            return(rep(0, ncol(yin_active_mat)))
          }
          grad_val <- frechet::SpheGeoGrad(yin_active_mat[i,], y_sphere)
          2 * dist_val * grad_val
        }))
        g_val <- colSums(grads_i_vals * s_in_weights_active)
        
        hess_terms_vals <- lapply(1:nrow(yin_active_mat), function(i) {
          dist_val <- frechet::SpheGeoDist(yin_active_mat[i,], y_sphere)
          if (abs(dist_val) < 1e-8 || abs(dist_val - pi) < 1e-8) { # Avoid issues at 0 or pi
            return(matrix(0, ncol(yin_active_mat), ncol(yin_active_mat)))
          }
          grad_val <- frechet::SpheGeoGrad(yin_active_mat[i,], y_sphere)
          hess_val <- frechet::SpheGeoHess(yin_active_mat[i,], y_sphere)
          term <- 2 * (grad_val %*% t(grad_val) + dist_val * hess_val)
          s_in_weights_active[i] * term
        })
        h_val <- Reduce("+", hess_terms_vals)
        return(list(value=f_val, gradient=g_val, hessian=h_val))
      }
      res_trust <- try(trust::trust(objFctn_sphere, y0, rinit=0.1, rmax=100, iterlim=100), silent=TRUE) # Increased rmax
      if (inherits(res_trust, "try-error") || !res_trust$converged) {
        warning("Trust algorithm failed for weighted spherical mean. Using initial guess.")
        l_hat_mean <- y0
      } else {
        l_hat_mean <- res_trust$argument / frechet:::l2norm(res_trust$argument)
      }
    } else { stop(paste("Unsupported metric space type:", metric_space_type)) }
  } # end normal case
  
  if (metric_space_type %in% c("covariance", "correlation")) { # Common post-processing for cov/corr
    l_hat_mat <- as.matrix(Matrix::nearPD(l_hat_mat_unproj, corr = (metric_space_type == "correlation"), keepDiag = (metric_space_type != "correlation"))$mat)
    l_hat_mat <- (l_hat_mat + t(l_hat_mat))/2
    if (metric_space_type == "correlation") {
      # Ensure diagonal is 1 after nearPD if it's a correlation matrix
      # nearPD with corr=TRUE should already do this, but enforce.
      # Also, ensure it remains positive semi-definite if possible.
      # A more careful projection to correlation matrices might be needed if nearPD doesn't perfectly achieve it.
      if(any(diag(l_hat_mat) <= 0)) { # Prevent sqrt of non-positive
        warning("Diagonal elements of estimated correlation matrix became non-positive after nearPD. Adjusting.")
        l_hat_mat <- cov2cor(Matrix::nearPD(l_hat_mat_unproj, corr=FALSE, keepDiag=FALSE, ensureSymmetry=TRUE, maxit=100)$mat) # Fallback: treat as cov, then convert
      } else {
        diag_inv_sqrt <- diag(1/sqrt(diag(l_hat_mat)), nrow=nrow(l_hat_mat))
        l_hat_mat <- diag_inv_sqrt %*% l_hat_mat %*% diag_inv_sqrt
        l_hat_mat <- (l_hat_mat + t(l_hat_mat))/2 # Re-symmetrize
      }
      diag(l_hat_mat) <- 1 # Force diagonal to be 1
      # Final check for PD, as scaling can break it.
      # If not PD, could try one more nearPD.
      if(min(eigen(l_hat_mat, symmetric=TRUE, only.values=TRUE)$values) < -1e-8) {
        l_hat_mat_corr <- try(as.matrix(Matrix::nearPD(l_hat_mat, corr=TRUE, ensureSymmetry=TRUE, keepDiag=TRUE, maxit=100)$mat), silent=TRUE)
        if(!inherits(l_hat_mat_corr, "try-error")) l_hat_mat <- l_hat_mat_corr
      }
    }
    l_hat_mean <- l_hat_mat
  } # else for density/sphere, l_hat_mean is already set
  
  return(l_hat_mean)
}


# Helper function to estimate one-sided Fréchet mean and variance
.estimate_one_sided_frechet_quantities <- function(Y_obj_list, X_scalar, metric_space_type,
                                                   c_val, h_frechet, K_frechet_fun, side,
                                                   frechet_optns, dist_fun_sq_calculator) {
  n_total <- length(Y_obj_list)
  s_in_weights <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_frechet, K_frechet_fun, side, n_total)
  
  active_indices <- which(s_in_weights != 0)
  if (length(active_indices) == 0) {
    warning(paste("No data points have non-zero weight for side", side, ". Cannot estimate Fréchet mean/variance."))
    return(list(l_hat = NULL, V_hat = NA, sigma_V_sq_hat = NA))
  }
  
  s_in_weights_active <- s_in_weights[active_indices]
  
  l_hat <- .estimate_weighted_frechet_mean(Y_obj_list, active_indices, s_in_weights_active, 
                                           metric_space_type, frechet_optns)
  if (is.null(l_hat)) {
    return(list(l_hat = NULL, V_hat = NA, sigma_V_sq_hat = NA)) # Error propagated
  }
  
  d_sq_vals_all_n <- sapply(1:n_total, function(i) {
    if (s_in_weights[i] == 0) return(0)
    dist_fun_sq_calculator(Y_obj_list[[i]], l_hat, frechet_optns)
  })
  
  V_hat <- sum(s_in_weights * d_sq_vals_all_n) / n_total
  
  term1_sigma_V <- sum(s_in_weights * (d_sq_vals_all_n^2)) / n_total # Sum s_in d^4 / n
  term2_sigma_V <- V_hat                                             # Sum s_in d^2 / n
  sigma_V_sq_hat <- term1_sigma_V - (term2_sigma_V^2)
  
  if (sigma_V_sq_hat < 0) {
    # warning(paste("Estimated sigma_V_sq_hat is negative for side", side, "(", sigma_V_sq_hat, "). Setting to small positive."))
    sigma_V_sq_hat <- max(sigma_V_sq_hat, 1e-12) # Allow very small negative due to precision, clamp harder if problem persists
  }
  if (sigma_V_sq_hat < 1e-12) sigma_V_sq_hat <- 1e-12 # Ensure positive for division
  
  return(list(l_hat = l_hat, V_hat = V_hat, sigma_V_sq_hat = sigma_V_sq_hat))
}

#' @title Test for Jumps in Metric-Space Conditional Means
#' @description Implements the ANOVA-style test for discontinuity in conditional Fréchet means.
#'
#' @param Y_obj A list of observed metric objects. Length n.
#' @param X_scalar A numeric vector of scalar covariates. Length n.
#' @param c_val The cutoff point for the discontinuity test.
#' @param metric_space_type A character string specifying the type of metric space.
#'        Supported: "density", "covariance", "correlation", "sphere".
#' @param h_frechet Numeric bandwidth for local Fréchet regression.
#' @param kernel_frechet_char Character, kernel for local Fréchet regression (e.g., "gauss", "epan").
#' @param frechet_optns A list of options for Fréchet-specific computations:
#'        - For "density": `qSup` (vector, required probability grid for quantiles),
#'          `den_opts_for_create_density` (list, options for `frechet::CreateDensity`, e.g., `kernelDen`, `bwDen`),
#'          `lower`, `upper` (optional, bounds for quantile function projection).
#'        - For "covariance", "correlation": `metric` (e.g., "frobenius", "power"), `alpha` (if metric="power").
#'        - For "sphere": (none specific beyond local regression params).
#' @param h_fx Numeric bandwidth for estimating f_X(c) (density of X at c). If NULL, uses `stats::bw.nrd0`.
#' @param kernel_fx_char Character, kernel for f_X(c) estimation (e.g., "gauss").
#'
#' @return A list containing the test statistic `Tn`, p-value `p_value`, and intermediate quantities.
#' @export
#' @examples
#' \dontrun{
#' # --- Example for Density Space (Wasserstein) ---
#' set.seed(123)
#' n1 <- 50; n2 <- 50; n_total <- n1 + n2
#' X_scalar <- c(runif(n1, 0, 0.5-0.01), runif(n2, 0.5+0.01, 1))
#' Y_obj_density <- vector("list", n_total)
#' for(i in 1:n_total) {
#'   if (X_scalar[i] < 0.5) {
#'     Y_obj_density[[i]] <- rnorm(100, mean = 0, sd = 1)
#'   } else {
#'     Y_obj_density[[i]] <- rnorm(100, mean = 0.5, sd = 1) # Jump in mean
#'   }
#' }
#' frechet_options_density <- list(
#'   qSup = seq(0.01, 0.99, length.out = 101),
#'   den_opts_for_create_density = list(kernelDen = "gauss") # bwDen will be auto by CreateDensity
#' )
#' test_result_density <- FrechetJumpTest(
#'   Y_obj = Y_obj_density, X_scalar = X_scalar, c_val = 0.5,
#'   metric_space_type = "density",
#'   h_frechet = 0.1, kernel_frechet_char = "epan",
#'   frechet_optns = frechet_options_density,
#'   h_fx = NULL, kernel_fx_char = "gauss"
#' )
#' print(test_result_density)
#'
#' # --- Example for Covariance Space (Frobenius) ---
#' set.seed(456)
#' p_cov <- 5 # Dimension of covariance matrices
#' Y_obj_cov <- vector("list", n_total)
#' true_mean_cov1 <- diag(p_cov)
#' true_mean_cov2 <- diag(p_cov) * 1.5 # Jump in scale
#' for(i in 1:n_total) {
#'   base_cov <- if (X_scalar[i] < 0.5) true_mean_cov1 else true_mean_cov2
#'   # Simulate sample cov matrix around base_cov (simplified)
#'   sample_data <- MASS::mvrnorm(20, mu = rep(0, p_cov), Sigma = base_cov)
#'   Y_obj_cov[[i]] <- cov(sample_data)
#' }
#' frechet_options_cov <- list(metric = "frobenius")
#' test_result_cov <- FrechetJumpTest(
#'   Y_obj = Y_obj_cov, X_scalar = X_scalar, c_val = 0.5,
#'   metric_space_type = "covariance",
#'   h_frechet = 0.1, kernel_frechet_char = "epan",
#'   frechet_optns = frechet_options_cov
#' )
#' print(test_result_cov)
#' }
FrechetJumpTest <- function(Y_obj, X_scalar, c_val, metric_space_type, 
                            h_frechet, 
                            kernel_frechet_char = "gauss", 
                            frechet_optns = list(), 
                            h_fx = NULL, 
                            kernel_fx_char = "gauss"
) {
  
  n_total <- length(Y_obj)
  if (n_total != length(X_scalar)) stop("Y_obj and X_scalar must have the same length.")
  if (h_frechet <= 0) stop("h_frechet (bandwidth for local Frechet regression) must be positive.")
  
  K_frechet_fun <- frechet:::kerFctn(kernel_frechet_char)
  
  dist_fun_sq_calculator <- function(y1, y2, opts) {
    if (metric_space_type == "density") {
      qSup_common <- opts$qSup %||% stop("qSup missing in frechet_optns for density distance")
      den_opts_cfd <- opts$den_opts_for_create_density %||% list()
      
      q1_vals <- .get_quantiles_for_obs_jump_test(y1, qSup_target = qSup_common, den_opts = den_opts_cfd)
      q2_vals <- .get_quantiles_for_obs_jump_test(y2, qSup_target = qSup_common, den_opts = den_opts_cfd)
      return(pracma::trapz(qSup_common, (q1_vals - q2_vals)^2))
    } else if (metric_space_type %in% c("covariance", "correlation")) {
      # y1, y2 are matrices, opts contain metric, alpha
      return(frechet::dist4cov(A = y1, B = y2, optns = opts)$dist^2)
    } else if (metric_space_type == "sphere") {
      # y1, y2 are vectors
      return(frechet::SpheGeoDist(y1, y2)^2)
    } else {
      stop(paste("Unsupported metric_space_type for dist_fun_sq:", metric_space_type))
    }
  }
  
  res_plus <- .estimate_one_sided_frechet_quantities(Y_obj, X_scalar, metric_space_type, c_val, h_frechet, K_frechet_fun, "+", frechet_optns, dist_fun_sq_calculator)
  if (is.null(res_plus$l_hat)) return(list(Tn = NA, p_value = NA, error = "Right side estimation failed"))
  
  res_minus <- .estimate_one_sided_frechet_quantities(Y_obj, X_scalar, metric_space_type, c_val, h_frechet, K_frechet_fun, "-", frechet_optns, dist_fun_sq_calculator)
  if (is.null(res_minus$l_hat)) return(list(Tn = NA, p_value = NA, error = "Left side estimation failed"))
  
  s_in_plus_weights <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_frechet, K_frechet_fun, "+", n_total)
  s_in_minus_weights <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_frechet, K_frechet_fun, "-", n_total)
  s_in_pooled_weights <- 0.5 * (s_in_plus_weights + s_in_minus_weights)
  
  active_indices_pooled <- which(s_in_pooled_weights != 0)
  l_hat_pooled <- .estimate_weighted_frechet_mean(Y_obj, active_indices_pooled, s_in_pooled_weights[active_indices_pooled], metric_space_type, frechet_optns)
  if (is.null(l_hat_pooled)) return(list(Tn = NA, p_value = NA, error = "Pooled mean estimation failed"))
  
  d_sq_vals_pooled_all_n <- sapply(1:n_total, function(i) {
    if (s_in_pooled_weights[i] == 0) return(0)
    dist_fun_sq_calculator(Y_obj[[i]], l_hat_pooled, frechet_optns)
  })
  V_hat_pooled <- sum(s_in_pooled_weights * d_sq_vals_pooled_all_n) / n_total
  
  if (is.null(h_fx)) h_fx <- stats::bw.nrd0(X_scalar)
  if (h_fx <=0) {
    warning("Bandwidth for f_X(c) (h_fx) is non-positive. Using a small default.")
    h_fx <- 0.1 * stats::sd(X_scalar) * (n_total^(-1/5)) # Silverman's rule of thumb style if bw.nrd0 failed
    if (h_fx <= 0) h_fx <- 1e-3 # Absolute fallback
  }
  
  fx_den_obj <- try(stats::density(X_scalar, bw = h_fx, kernel = kernel_fx_char, n = 1024, from = min(X_scalar)-3*h_fx, to = max(X_scalar)+3*h_fx), silent=TRUE)
  if(inherits(fx_den_obj, "try-error")){
    warning("stats::density for f_X(c) failed. Using uniform as fallback for f_X_hat_c.")
    f_X_hat_c <- 1 / (max(X_scalar) - min(X_scalar))
    if(!is.finite(f_X_hat_c) || f_X_hat_c <=0) f_X_hat_c <- 1
  } else {
    f_X_hat_c <- stats::approx(fx_den_obj$x, fx_den_obj$y, xout = c_val)$y
  }
  
  if (is.na(f_X_hat_c) || f_X_hat_c < 1e-9) {
    warning(paste("Estimated f_X(c) is problematic (value:", f_X_hat_c, "). Clamping to 1e-9."))
    f_X_hat_c <- 1e-9
  }
  
  integrate_K_moment_internal <- function(K_fun_char_internal, k_pow, j_pow) {
    K_to_integrate <- function(u_vec) {
      sapply(u_vec, function(u) {
        if (u < 0) return(0)
        val_K <- frechet:::kerFctn(K_fun_char_internal)(u)
        u^j_pow * val_K^k_pow
      })
    }
    limit_val <- if (K_fun_char_internal %in% c("rect", "epan", "quar")) 1 else 8
    # Catch cases where integration range is empty or problematic for specific kernels
    if (limit_val == 0 && K_fun_char_internal != "gauss") { # e.g. for rect kernel, integral from 0 to 0 if limit would be 0
      if (K_fun_char_internal == "rect" && j_pow==0 && k_pow==1) return(0.5) # K(0)*0 - (-K(0)*0) half of integral K(u)du from 0 to 1 for rect
      # this needs more careful thought for specific K moments at u=0 if limit is 0
    }
    
    res <- try(stats::integrate(K_to_integrate, lower = 0, upper = limit_val, subdivisions = 2000, rel.tol = 1e-7, stop.on.error = FALSE)$value, silent=TRUE)
    if (inherits(res, "try-error") || is.na(res)) {
      warning(paste("Integration failed for K_moment:", K_fun_char_internal, k_pow, j_pow, ". Result:", res))
      return(NA)
    }
    return(res)
  }
  
  K_p_10 <- integrate_K_moment_internal(kernel_frechet_char, 1, 0)
  K_p_11 <- integrate_K_moment_internal(kernel_frechet_char, 1, 1)
  K_p_12 <- integrate_K_moment_internal(kernel_frechet_char, 1, 2)
  
  integrand_S_num_internal <- function(u_vec) {
    sapply(u_vec, function(u) {
      if (u < 0) return(0)
      K_val <- frechet:::kerFctn(kernel_frechet_char)(u)
      term_in_paren <- K_p_12 - u * K_p_11
      (term_in_paren^2) * (K_val^2)
    })
  }
  limit_S_num_internal <- if (kernel_frechet_char %in% c("rect", "epan", "quar")) 1 else 8
  S_num_val <- try(stats::integrate(integrand_S_num_internal, lower = 0, upper = limit_S_num_internal, subdivisions = 2000, rel.tol = 1e-7, stop.on.error = FALSE)$value, silent=TRUE)
  S_den_val_paren <- (K_p_12 * K_p_10 - K_p_11^2)
  
  if (inherits(S_num_val, "try-error") || is.na(S_num_val) || abs(S_den_val_paren) < 1e-9) {
    warning("Could not compute S_K constant. Check kernel and its moments.")
    S_K <- NA
  } else {
    S_K <- S_num_val / (S_den_val_paren^2)
  }
  
  if(is.na(S_K) || !is.finite(S_K)) {
    warning("S_K is NA or Inf. Test statistic will be NA.")
    return(list(Tn=NA, p_value=NA, error="S_K computation failed"))
  }
  
  F_n <- V_hat_pooled - 0.5 * (res_plus$V_hat + res_minus$V_hat)
  
  U_n_numerator <- (res_plus$V_hat - res_minus$V_hat)^2
  U_n_denominator <- (S_K * res_plus$sigma_V_sq_hat / f_X_hat_c) + (S_K * res_minus$sigma_V_sq_hat / f_X_hat_c)
  
  U_n <- if (abs(U_n_denominator) < 1e-12) {
    if (abs(U_n_numerator) < 1e-12) 0 else Inf 
  } else {
    U_n_numerator / U_n_denominator
  }
  
  T_n_term1 <- n_total * h_frechet * U_n
  T_n_term2 <- if (abs(U_n_denominator) < 1e-12) {
    if (abs(F_n^2) < 1e-12) 0 else Inf
  } else {
    (n_total * h_frechet * F_n^2) / U_n_denominator
  }
  
  T_n <- T_n_term1 + T_n_term2
  
  p_value <- if (is.na(T_n) || is.infinite(T_n) || T_n < 0) { # T_n should be non-negative
    NA 
  } else {
    stats::pchisq(T_n, df = 1, lower.tail = FALSE)
  }
  
  return(list(Tn = T_n, p_value = p_value,
              V_hat_plus = res_plus$V_hat, V_hat_minus = res_minus$V_hat, V_hat_pooled = V_hat_pooled,
              sigma_V_sq_hat_plus = res_plus$sigma_V_sq_hat, sigma_V_sq_hat_minus = res_minus$sigma_V_sq_hat,
              f_X_hat_c = f_X_hat_c, S_K = S_K,
              F_n = F_n, U_n = U_n,
              l_hat_plus = res_plus$l_hat, l_hat_minus = res_minus$l_hat, l_hat_pooled = l_hat_pooled,
              kernel_constants = list(K_p_10=K_p_10, K_p_11=K_p_11, K_p_12=K_p_12)))
}