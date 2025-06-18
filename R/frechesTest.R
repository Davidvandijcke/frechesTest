# R/frechesTest.R

#' Test for Jumps in Metric-Space Conditional Frechet Means
#'
#' Implements an ANOVA-style test for detecting discontinuities (jumps) in the 
#' conditional Frechet mean of random objects valued in arbitrary metric spaces.
#' The test is based on the methodology described in "A Test for Jumps in 
#' Metric-Space Conditional Means" by David Van Dijcke.
#'
#' @param Y_obj A list of observed metric objects. Each element corresponds to one
#'   observation. For \code{metric_space_type = "density"}, elements can be raw 
#'   numeric vectors, histogram objects, or specific quantile lists. For 
#'   \code{"covariance"} or \code{"correlation"}, elements should be matrices. 
#'   For \code{"sphere"}, elements should be numeric unit vectors. For 
#'   \code{"network"}, elements should be graph Laplacian matrices.
#' @param X_scalar A numeric vector of scalar covariates, corresponding to \code{Y_obj}.
#'   Must have the same length as \code{Y_obj}.
#' @param c_val The scalar cutoff point in \code{X_scalar} where the discontinuity 
#'   is tested, splitting the data into \code{X_scalar >= c} and \code{X_scalar < c}.
#' @param metric_space_type A character string specifying the type of metric space.
#'   Supported types: \code{"density"} (for probability distributions using 
#'   Wasserstein-2 distance), \code{"covariance"} (for SPD matrices), 
#'   \code{"correlation"} (for correlation matrices), \code{"sphere"} (for unit 
#'   vectors on hypersphere), \code{"network"} (for graph Laplacians).
#' @param h_frechet Bandwidth for local Frechet regression. Can be:
#'   \itemize{
#'     \item Numeric scalar: Fixed bandwidth value
#'     \item \code{"CV"}: Use K-fold cross-validation for automatic selection
#'   }
#' @param kernel_frechet_char Character string specifying the kernel function for 
#'   local Frechet regression weights. Options include \code{"gauss"}, \code{"epan"}, 
#'   \code{"rect"}, \code{"quar"}. Default is \code{"epan"}.
#' @param frechet_optns A list of options specific to the \code{metric_space_type}:
#'   \itemize{
#'     \item For \code{"density"}:
#'       \itemize{
#'         \item \code{qSup}: (Required) Numeric vector spanning [0,1] for quantile grid
#'         \item \code{den_opts_for_create_density}: Options for \code{frechet::CreateDensity}
#'         \item \code{lower}, \code{upper}: Bounds for quantile function support
#'         \item \code{osqp_model_for_density_proj}: Cached OSQP model for optimization
#'       }
#'     \item For \code{"covariance"} or \code{"correlation"}:
#'       \itemize{
#'         \item \code{metric}: Distance metric (\code{"frobenius"}, \code{"power"}, 
#'               \code{"log_cholesky"}, \code{"cholesky"})
#'         \item \code{alpha}: Power parameter for \code{"power"} metric
#'       }
#'     \item For \code{"network"}:
#'       \itemize{
#'         \item \code{network_directed}: Logical, whether graphs are directed. Default FALSE.
#'         \item \code{W_laplacian_bound}: Upper bound for edge weights
#'         \item \code{osqp_model_for_laplacian_proj}: Cached OSQP model
#'       }
#'      \item For \code{"sphere"}:
#'      \itemize{
#'      \item \code{enforce_positive}: Logical, whether to force Frechet mean into
#'      positive orthant. Default FALSE. 
#'      }
#'      }
#' @param h_fx Numeric scalar, bandwidth for estimating the density of \code{X_scalar} 
#'   at \code{c_val}. If \code{NULL}, uses \code{stats::bw.nrd0(X_scalar)}.
#' @param kernel_fx_char Character string, kernel for density estimation of X. 
#'   Default is \code{"gauss"}.
#' @param cv_K_folds Integer, number of folds for cross-validation when 
#'   \code{h_frechet = "CV"}. Default is 5.
#' @param cv_n_bw_candidates Integer, number of bandwidth candidates to test 
#'   during cross-validation. Default is 10.
#' @param undersmooth_factor Numeric, factor for undersmoothing the selected 
#'   bandwidth for variance estimation. Default is 0.9.
#' @param min_eff_points_cv_fold Integer, minimum effective sample size required 
#'   in each CV fold. Default is 5.
#' @param min_bw_cv_factor Numeric, minimum bandwidth as fraction of data range. 
#'   Default is 0.001.
#' @param verbose Logical, whether to print progress messages. Default is \code{FALSE}.
#'
#' @return A list containing:
#'   \item{Tn}{The test statistic value}
#'   \item{p_value}{P-value based on Chi-squared distribution with 1 df}
#'   \item{V_hat_plus}{Estimated conditional Frechet variance from the right of \code{c_val}}
#'   \item{V_hat_minus}{Estimated conditional Frechet variance from the left of \code{c_val}}
#'   \item{V_hat_pooled}{Estimated pooled conditional Frechet variance around \code{c_val}}
#'   \item{sigma_V_sq_hat_plus}{Estimated variance of squared distances from the right}
#'   \item{sigma_V_sq_hat_minus}{Estimated variance of squared distances from the left}
#'   \item{f_X_hat_c}{Estimated density of X at cutoff point}
#'   \item{S_K}{Kernel constant for asymptotic distribution}
#'   \item{lambda_hat}{Estimated mixing proportion}
#'   \item{F_n}{Mean difference component of test statistic}
#'   \item{U_n}{Variance difference component of test statistic}
#'   \item{l_hat_plus}{Estimated Frechet mean from the right}
#'   \item{l_hat_minus}{Estimated Frechet mean from the left}
#'   \item{l_hat_pooled}{Estimated pooled Frechet mean}
#'   \item{kernel_constants}{List of computed kernel moments}
#'   \item{num_active_plus}{Number of observations with non-zero weights (right side)}
#'   \item{num_active_minus}{Number of observations with non-zero weights (left side)}
#'   \item{num_active_pooled}{Number of observations with non-zero pooled weights}
#'   \item{h_mean_cv_selected}{Selected bandwidth (if CV was used)}
#'   \item{h_variance_used}{Bandwidth used for variance estimation}
#'   \item{D_n_U_denominator}{Denominator of U component}
#'   \item{D_n_F_denominator}{Denominator of F component}
#'   \item{error}{Error message if computation failed, otherwise NULL}
#'
#' @details
#' The test combines differences in estimated Frechet means and variances from both 
#' sides of the cutoff. The test statistic follows a Chi-squared distribution with 
#' 1 degree of freedom under the null hypothesis of no jump.
#'
#' For \code{metric_space_type = "density"}, observations can be:
#' \itemize{
#'   \item Numeric vectors (raw data)
#'   \item Lists with \code{q} (quantile values), \code{qSup} (support), and \code{type="quantile"}
#'   \item Histogram objects with \code{mids}, \code{counts}, \code{breaks}
#' }
#'
#' The \code{frechet_optns$qSup} parameter is essential for density space computations
#' as it defines the common grid for L2-Wasserstein distance calculations.
#'
#' @examples
#' \dontrun{
#' # Example 1: Density space with jump in mean
#' set.seed(123)
#' n <- 100
#' X <- runif(n)
#' Y <- vector("list", n)
#' 
#' for(i in 1:n) {
#'   if (X[i] < 0.5) {
#'     Y[[i]] <- rnorm(50, mean = 0, sd = 1)
#'   } else {
#'     Y[[i]] <- rnorm(50, mean = 2, sd = 1)  # Jump in mean
#'   }
#' }
#' 
#' result_density <- frechesTest(
#'   Y_obj = Y,
#'   X_scalar = X,
#'   c_val = 0.5,
#'   metric_space_type = "density",
#'   h_frechet = "CV",
#'   frechet_optns = list(qSup = seq(0.01, 0.99, length.out = 99))
#' )
#' 
#' print(paste("P-value:", result_density$p_value))
#' 
#' # Example 2: Covariance space
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   Y_cov <- vector("list", n)
#'   for(i in 1:n) {
#'     if (X[i] < 0.5) {
#'       Sigma <- diag(3)
#'     } else {
#'       Sigma <- diag(c(2, 2, 2))  # Jump in scale
#'     }
#'     data <- MASS::mvrnorm(30, mu = rep(0, 3), Sigma = Sigma)
#'     Y_cov[[i]] <- cov(data)
#'   }
#'   
#'   result_cov <- frechesTest(
#'     Y_obj = Y_cov,
#'     X_scalar = X,
#'     c_val = 0.5,
#'     metric_space_type = "covariance",
#'     h_frechet = 0.2,
#'     frechet_optns = list(metric = "frobenius")
#'   )
#'   
#'   print(paste("P-value:", result_cov$p_value))
#' }
#' }
#'
#' @references
#' Van Dijcke, D. (2023). "A Test for Jumps in Metric-Space Conditional Means." 
#' Working Paper, University of Michigan.
#'
#' @importFrom Matrix Matrix nearPD sparseMatrix
#' @importFrom fdadensity dens2quantile
#' @importFrom frechet CreateDensity SpheGeoDist SpheGeoGrad SpheGeoHess dist4cov 
#' @importFrom osqp osqp osqpSettings
#' @importFrom pracma trapz
#' @importFrom stats IQR approx bw.nrd0 cov density integrate pchisq sd
#' @importFrom trust trust
#' @export
frechesTest <- function(Y_obj, X_scalar, c_val, metric_space_type,
                        h_frechet = "CV",
                        kernel_frechet_char = "epan",
                        frechet_optns = list(),
                        h_fx = NULL,
                        kernel_fx_char = "gauss",
                        cv_K_folds = 5,
                        cv_n_bw_candidates = 10,
                        undersmooth_factor = 0.9,
                        min_eff_points_cv_fold = 5, # Min effective points for Frechet mean in CV fold
                        min_bw_cv_factor = 0.001,
                        verbose = FALSE
) {
  
  # --- Input Validations & Setup ---
  
  # Validate basic inputs
  if (!is.list(Y_obj)) {
    stop("Y_obj must be a list of metric objects.")
  }
  
  if (!is.numeric(X_scalar)) {
    stop("X_scalar must be a numeric vector.")
  }
  
  if (length(Y_obj) != length(X_scalar)) {
    stop("Y_obj and X_scalar must have the same length. Y_obj has ", 
         length(Y_obj), " elements, X_scalar has ", length(X_scalar), " elements.")
  }
  if (metric_space_type == "sphere" && isTRUE(frechet_optns$enforce_positive)) {
    any_neg <- any(sapply(Y_obj, function(v) any(v < 0)))
    if (any_neg)
      warning("enforce_positive = TRUE but some input vectors have negative coordinates.")
  }
  
  
  if (length(Y_obj) == 0) {
    stop("Y_obj and X_scalar cannot be empty.")
  }
  
  if (!is.numeric(c_val) || length(c_val) != 1) {
    stop("c_val must be a single numeric value.")
  }
  
  if (!is.finite(c_val)) {
    stop("c_val must be finite.")
  }
  
  # Validate metric_space_type
  valid_spaces <- c("density", "covariance", "correlation", "sphere", "network")
  if (!is.character(metric_space_type) || length(metric_space_type) != 1 || 
      !metric_space_type %in% valid_spaces) {
    stop("metric_space_type must be one of: ", paste(valid_spaces, collapse = ", "))
  }
  
  # Validate h_frechet
  if (is.character(h_frechet)) {
    if (length(h_frechet) != 1 || toupper(h_frechet) != "CV") {
      stop("If h_frechet is character, it must be 'CV' for cross-validation.")
    }
  } else if (is.numeric(h_frechet)) {
    if (length(h_frechet) != 1 || h_frechet <= 0 || !is.finite(h_frechet)) {
      stop("If h_frechet is numeric, it must be a single positive finite value.")
    }
  } else {
    stop("h_frechet must be either 'CV' or a positive numeric value.")
  }
  
  # Validate CV parameters if using CV
  if (is.character(h_frechet) && toupper(h_frechet) == "CV") {
    if (!is.numeric(cv_K_folds) || length(cv_K_folds) != 1 || 
        cv_K_folds < 2 || cv_K_folds != round(cv_K_folds)) {
      stop("cv_K_folds must be an integer >= 2.")
    }
    
    if (!is.numeric(cv_n_bw_candidates) || length(cv_n_bw_candidates) != 1 || 
        cv_n_bw_candidates < 2 || cv_n_bw_candidates != round(cv_n_bw_candidates)) {
      stop("cv_n_bw_candidates must be an integer >= 2.")
    }
    
    if (!is.numeric(undersmooth_factor) || length(undersmooth_factor) != 1 || 
        undersmooth_factor <= 0 || undersmooth_factor > 1) {
      stop("undersmooth_factor must be a numeric value in (0, 1].")
    }
  }
  
  # Validate kernel
  if (!is.character(kernel_frechet_char) || length(kernel_frechet_char) != 1) {
    stop("kernel_frechet_char must be a single character string.")
  }
  
  # Validate h_fx if provided
  if (!is.null(h_fx)) {
    if (!is.numeric(h_fx) || length(h_fx) != 1 || h_fx <= 0 || !is.finite(h_fx)) {
      stop("h_fx must be NULL or a single positive finite numeric value.")
    }
  }
  
  # Validate frechet_optns
  if (!is.list(frechet_optns)) {
    stop("frechet_optns must be a list.")
  }
  
  # Check if we have enough data on both sides of cutoff
  n_left <- sum(X_scalar < c_val)
  n_right <- sum(X_scalar >= c_val)
  
  if (n_left == 0) {
    stop("No observations found with X_scalar < c_val (", c_val, "). ",
         "Cannot perform jump test.")
  }
  
  if (n_right == 0) {
    stop("No observations found with X_scalar >= c_val (", c_val, "). ",
         "Cannot perform jump test.")
  }
  
  min_side_size <- 3  # Minimum observations per side for meaningful estimation
  if (n_left < min_side_size) {
    warning("Only ", n_left, " observations with X_scalar < c_val. ",
            "Results may be unreliable with fewer than ", min_side_size, " observations per side.")
  }
  
  if (n_right < min_side_size) {
    warning("Only ", n_right, " observations with X_scalar >= c_val. ",
            "Results may be unreliable with fewer than ", min_side_size, " observations per side.")
  }
  if (metric_space_type == "density") {
    if (is.null(frechet_optns$qSup)) {
      stop("For metric_space_type='density', frechet_optns$qSup must be provided.")
    }
    qSup_val <- frechet_optns$qSup
    # Basic validation for qSup
    if (!is.numeric(qSup_val) || any(is.na(qSup_val)) || 
        abs(min(qSup_val) - 0) > 1e-9 || abs(max(qSup_val) - 1) > 1e-9) { # Allow small tolerance for 0 and 1
      stop("For metric_space_type='density', frechet_optns$qSup must be a numeric vector spanning [0, 1].")
    }
    if(any(diff(qSup_val) <= 0)) {
      stop("For metric_space_type='density', frechet_optns$qSup must be strictly increasing.")
    }
    
    # --- Default for network_directed and ensure it's in frechet_optns ---
    is_network_directed <- frechet_optns$network_directed %||% FALSE
    if (metric_space_type == "network" && is.null(frechet_optns$network_directed)) {
      if (verbose) {
        cat("Note: frechet_optns$network_directed not specified, defaulting to FALSE (undirected graph Laplacians).\n")
      }
    }
    frechet_optns$network_directed <- is_network_directed # Ensure flag is consistently available
    
    ## OPTIMIZATION: Pre-calculate quantiles for all density objects
    qSup_reference <- frechet_optns$qSup
    den_opts_reference <- frechet_optns$den_opts_for_create_density %||% list()
    
    if (verbose) cat("Pre-calculating quantiles for density objects...\n")
    
    Y_obj_processed <- vector("list", length(Y_obj))
    any_errors_precomputation <- FALSE
    
    is_precomputed_quantile <- function(obj, target_qSup) {
      is.list(obj) && !is.null(obj$type) && obj$type == "quantile" &&
        !is.null(obj$q) && !is.null(obj$qSup) &&
        length(obj$qSup) == length(target_qSup) && all(abs(obj$qSup - target_qSup) < 1e-9)
    }
    
    for (i in seq_along(Y_obj)) {
      if (is_precomputed_quantile(Y_obj[[i]], qSup_reference)) {
        Y_obj_processed[[i]] <- Y_obj[[i]]
      } else {
        quantiles_vec_i <- tryCatch({
          .get_quantiles_for_obs_jump_test(obs = Y_obj[[i]],
                                           qSup_target = qSup_reference,
                                           den_opts = den_opts_reference)
        }, error = function(e) {
          warning(paste0("Error pre-calculating quantiles for Y_obj[[", i, "]]: ", e$message, 
                         ". Original object was of class: ", class(Y_obj[[i]])[1]))
          any_errors_precomputation <<- TRUE # Use <<- to assign to parent env var
          NULL
        })
        
        if (is.null(quantiles_vec_i)) {
          Y_obj_processed[[i]] <- NULL # Mark failure
        } else {
          Y_obj_processed[[i]] <- list(q = quantiles_vec_i, qSup = qSup_reference, type = "quantile")
        }
      }
    }
    
    if (any_errors_precomputation || any(sapply(Y_obj_processed, is.null))) {
      stop("Failed to pre-compute quantiles for one or more density objects. Check warnings and ensure frechet_optns$qSup is appropriate and data format is supported.")
    }
    Y_obj <- Y_obj_processed # Replace Y_obj with the list of pre-computed quantile objects
    if (verbose) cat("Pre-calculation of quantiles complete.\n")
    
    # Initialize cache field for density OSQP model if not present in frechet_optns
    if (is.null(frechet_optns$osqp_model_for_density_proj)) {
      frechet_optns$osqp_model_for_density_proj <- NULL 
    }
  } # End density pre-computation
  
  n_total <- length(Y_obj)
  if (n_total != length(X_scalar)) stop("Y_obj and X_scalar must have the same length.")
  
  # Ensure frechet::kerFctn is available
  if (!exists("kerFctn", mode = "function", envir = asNamespace("frechet"))) {
    stop("kerFctn from 'frechet' package components not found. Please ensure 'frechet' is loaded or its relevant components are sourced.")
  }
  actual_kernel_function_to_pass <- get("kerFctn", envir = asNamespace("frechet"))(kernel_frechet_char)
  
  # --- OSQP Model Setup for Laplacian Projection (if applicable) ---
  # This sets up frechet_optns$osqp_model_for_laplacian_proj which is then *used* by
  # .project_to_laplacian_space_osqp via .estimate_weighted_frechet_mean
  if (metric_space_type == "network") {
    # Check if model needs to be (re-)built: not cached, or M_nodes has changed
    M_nodes_proj <- ncol(Y_obj[[1]]) # Assuming Y_obj is not empty and Y_obj[[1]] is a matrix
    if(is.null(M_nodes_proj) || M_nodes_proj < 1){
      stop("Could not determine a valid number of nodes (M_nodes_proj) from Y_obj[[1]] for network projection.")
    }
    
    rebuild_network_osqp_model <- TRUE # Default to rebuild
    if (!is.null(frechet_optns$osqp_model_for_laplacian_proj) && 
        !is.null(frechet_optns$osqp_model_for_laplacian_proj$M_nodes) &&
        frechet_optns$osqp_model_for_laplacian_proj$M_nodes == M_nodes_proj) {
      rebuild_network_osqp_model <- FALSE # Cached model matches current M_nodes
    }
    
    if (rebuild_network_osqp_model) {
      if(verbose) cat("Building OSQP model for Laplacian projection (M_nodes =", M_nodes_proj, ")...\n")
      W_lap_bound_proj <- frechet_optns$W_laplacian_bound %||% 1.0
      
      P_osqp_sparse_proj <- Matrix::Matrix(diag(M_nodes_proj^2), sparse = TRUE)
      
      num_off_diag_pairs_proj <- M_nodes_proj * (M_nodes_proj - 1) / 2
      # Constraints: Symmetry, Zero row sums, Off-diagonal bounds
      # Total constraints: num_off_diag_pairs (sym) + M_nodes (rowsum) + num_off_diag_pairs (offdiag_bounds)
      num_constraints_sym_proj <- if(M_nodes_proj > 1) num_off_diag_pairs_proj else 0
      num_constraints_rowsum_proj <- M_nodes_proj
      num_constraints_offdiag_bounds_proj <- if(M_nodes_proj > 1) num_off_diag_pairs_proj else 0
      num_total_constraints_proj <- num_constraints_sym_proj + num_constraints_rowsum_proj + num_constraints_offdiag_bounds_proj
      
      if (num_total_constraints_proj == 0 && M_nodes_proj == 1) { # Special case: 1 node graph
        const_matrix_A_sparse_proj <- Matrix::Matrix(0, nrow = 1, ncol = 1, sparse = TRUE) # L11 = 0
        const_matrix_A_sparse_proj[1,1] <- 1
        l_bounds_proj <- 0
        u_bounds_proj <- 0
      } else if (num_total_constraints_proj > 0) {
        const_matrix_A_proj <- matrix(0, nrow = num_total_constraints_proj, ncol = M_nodes_proj^2)
        l_bounds_proj <- numeric(num_total_constraints_proj)
        u_bounds_proj <- numeric(num_total_constraints_proj)
        
        current_row_proj <- 0
        if (num_constraints_sym_proj > 0) {
          for (c_node in 2:M_nodes_proj) {
            for (r_node in 1:(c_node - 1)) {
              current_row_proj <- current_row_proj + 1
              idx_rc_in_vec <- (c_node - 1) * M_nodes_proj + r_node
              idx_cr_in_vec <- (r_node - 1) * M_nodes_proj + c_node
              const_matrix_A_proj[current_row_proj, idx_rc_in_vec] <- 1
              const_matrix_A_proj[current_row_proj, idx_cr_in_vec] <- -1
              l_bounds_proj[current_row_proj] <- 0; u_bounds_proj[current_row_proj] <- 0
            }
          }
        }
        
        for (r_node in 1:M_nodes_proj) { # Zero row sums
          current_row_proj <- current_row_proj + 1
          for (c_node in 1:M_nodes_proj) {
            idx_rc_in_vec <- (c_node - 1) * M_nodes_proj + r_node
            const_matrix_A_proj[current_row_proj, idx_rc_in_vec] <- 1
          }
          l_bounds_proj[current_row_proj] <- 0; u_bounds_proj[current_row_proj] <- 0
        }
        
        if (num_constraints_offdiag_bounds_proj > 0) {
          for (c_node in 2:M_nodes_proj) { # Off-diagonal bounds -W <= L_ij <= 0 for i < j
            for (r_node in 1:(c_node - 1)) {
              current_row_proj <- current_row_proj + 1
              idx_rc_in_vec <- (c_node - 1) * M_nodes_proj + r_node
              const_matrix_A_proj[current_row_proj, idx_rc_in_vec] <- 1
              l_bounds_proj[current_row_proj] <- -W_lap_bound_proj
              u_bounds_proj[current_row_proj] <- 0
            }
          }
        }
        const_matrix_A_sparse_proj <- Matrix::Matrix(const_matrix_A_proj, sparse = TRUE)
      } else { # Should not happen for M_nodes_proj > 1
        stop("Error in setting up constraints for network OSQP model.")
      }
      
      default_osqp_settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 20000, polish = TRUE)
      current_osqp_settings_proj <- frechet_optns$osqp_settings_laplacian_proj %||% default_osqp_settings
      
      # Use a placeholder q; it will be updated in .project_to_laplacian_space_osqp
      q_placeholder_network <- rep(0, M_nodes_proj^2) 
      osqp_model_obj_for_proj <- osqp::osqp(P = P_osqp_sparse_proj, q = q_placeholder_network,
                                            A = const_matrix_A_sparse_proj,
                                            l = l_bounds_proj, u = u_bounds_proj,
                                            pars = current_osqp_settings_proj)
      # Store the model and M_nodes in frechet_optns for reuse
      frechet_optns$osqp_model_for_laplacian_proj <- list(model = osqp_model_obj_for_proj, M_nodes = M_nodes_proj)
      if(verbose) cat("OSQP model for Laplacian projection built and cached.\n")
    } else {
      if(verbose) cat("Using cached OSQP model for Laplacian projection.\n")
    }
  }
  
  
  # --- dist_fun_sq_calculator (defined once) ---
  dist_fun_sq_calculator <- function(y1, y2, opts) { # opts is frechet_optns
    if (metric_space_type == "density") {
      qSup_common <- opts$qSup %||% stop("qSup common grid missing for density distance.")
      den_opts_cfd <- opts$den_opts_for_create_density %||% list()
      
      # y1 comes from Y_obj, so it should be a pre-computed quantile list
      # .get_quantiles_for_obs_jump_test will efficiently extract q1_vals
      q1_vals <- .get_quantiles_for_obs_jump_test(y1, qSup_target = qSup_common, den_opts = den_opts_cfd)
      
      # y2 is the Frechet mean, which is already a list(q=..., qSup=..., type="quantile")
      if (!is.list(y2) || is.null(y2$q) || !identical(y2$type, "quantile")) {
        stop("Internal error: y2 (Frechet mean) in dist_fun_sq_calculator (density) is not expected list structure {q, qSup, type='quantile'}.")
      }
      q2_numeric_vector <- y2$q # Extract the quantile vector
      
      if (length(q1_vals) != length(qSup_common) || length(q2_numeric_vector) != length(qSup_common)) {
        stop(paste("Internal error: Quantile vector lengths (q1:", length(q1_vals), ", q2:", length(q2_numeric_vector),
                   ") do not match qSup_common length (", length(qSup_common), ") in dist_fun_sq_calculator.", sep=""))
      }
      return(pracma::trapz(qSup_common, (q1_vals - q2_numeric_vector)^2))
    } else if (metric_space_type %in% c("covariance", "correlation", "network")) {
      # For network, y2 (Frechet mean) is a projected Laplacian matrix. y1 is an observed matrix.
      # frechet::dist4cov should handle this.
      return(frechet::dist4cov(A = y1, B = y2, optns = opts)$dist^2)
    } else if (metric_space_type == "sphere") {
      return(frechet::SpheGeoDist(y1, y2)^2)
    } else {
      stop(paste("Unsupported metric_space_type for dist_fun_sq:", metric_space_type))
    }
  }
  
  # --- Bandwidth Selection and Undersmoothing ---
  h_mean_cv_selected <- NA_real_
  h_variance_to_use <- NA_real_
  
  if (is.character(h_frechet) && toupper(h_frechet) == "CV") {
    if (verbose) cat("Performing K-Fold Cross-Validation for a common Frechet mean bandwidth for the jump...\n")
    
    h_mean_cv_selected <- ._cv_for_frechet_jump_h(
      Y_obj_list = Y_obj, X_scalar = X_scalar, c_val = c_val,
      metric_space_type = metric_space_type, kernel_frechet_char = kernel_frechet_char,
      actual_kernel_function = actual_kernel_function_to_pass,
      frechet_optns = frechet_optns, dist_fun_sq_calculator = dist_fun_sq_calculator,
      K_folds = cv_K_folds, n_bw_candidates = cv_n_bw_candidates,
      min_eff_points_cv_fold = min_eff_points_cv_fold,
      min_bw_heuristic_factor = min_bw_cv_factor
    )
    if (verbose) cat(paste("  CV selected common h_mean:", signif(h_mean_cv_selected,4), "\n"))
    
    h_variance_to_use <- h_mean_cv_selected * undersmooth_factor
    if (verbose) cat(paste("  Using h_variance (undersmoothed):", signif(h_variance_to_use,4), "\n"))
    
  } else if (is.numeric(h_frechet)) {
    if (length(h_frechet) != 1 || h_frechet <= 0) stop("User-supplied h_frechet must be a single positive number.")
    h_mean_cv_selected <- h_frechet # This is h_F in paper, used for means in CV
    h_variance_to_use <- h_mean_cv_selected * undersmooth_factor # This is h_V in paper, used for variances
    if (verbose) cat(paste("  User-supplied h_mean:", signif(h_mean_cv_selected,4), ", using h_variance (undersmoothed):", signif(h_variance_to_use,4), "\n"))
  } else {
    stop("h_frechet must be numeric or 'CV'.")
  }
  
  if (is.na(h_variance_to_use) || h_variance_to_use <= 1e-5) { # Increased robustness for very small/NA h_variance
    old_h_var <- h_variance_to_use
    # Use h_mean_cv_selected (which could be user-supplied or CV-result) as a basis
    h_mean_finite_base <- if(is.finite(h_mean_cv_selected) && h_mean_cv_selected > 1e-5) h_mean_cv_selected else {
      # Fallback for h_mean_cv_selected if it's bad
      sd_X_fallback <- stats::sd(X_scalar, na.rm=TRUE); if(is.na(sd_X_fallback) || sd_X_fallback < 1e-5) sd_X_fallback <- 1
      0.1 * sd_X_fallback
    }
    h_variance_to_use <- max(1e-5, h_mean_finite_base * 0.1, na.rm = TRUE) # Ensure it's positive and not excessively small
    warning(paste("Calculated h_variance was NA, non-positive or too small (", signif(old_h_var,4), "). Reset to ", signif(h_variance_to_use,4), 
                  " based on h_mean_cv_selected or fallback."))
  }
  
  # --- One-sided estimates using h_variance_to_use ---
  if (verbose) cat("Estimating one-sided Frechet quantities...\n")
  res_plus <- .estimate_one_sided_frechet_quantities(Y_obj, X_scalar, metric_space_type, c_val,
                                                     h_variance_to_use, # Use h_V for final variance estimates
                                                     actual_kernel_function_to_pass,
                                                     "+", frechet_optns, dist_fun_sq_calculator)
  if (is.null(res_plus$l_hat)) {
    return(list(Tn = NA, p_value = NA, error = "Right side Frechet estimation failed.",
                h_mean_cv_selected = h_mean_cv_selected, h_variance_used = h_variance_to_use,
                debug_info_plus = res_plus$debug_info, debug_info_minus = NULL))
  }
  
  res_minus <- .estimate_one_sided_frechet_quantities(Y_obj, X_scalar, metric_space_type, c_val,
                                                      h_variance_to_use, # Use h_V
                                                      actual_kernel_function_to_pass,
                                                      "-", frechet_optns, dist_fun_sq_calculator)
  if (is.null(res_minus$l_hat)) {
    return(list(Tn = NA, p_value = NA, error = "Left side Frechet estimation failed.",
                h_mean_cv_selected = h_mean_cv_selected, h_variance_used = h_variance_to_use,
                debug_info_plus = res_plus$debug_info, debug_info_minus = res_minus$debug_info))
  }
  
  # --- Pooled estimates using h_variance_to_use ---
  if (verbose) cat("Estimating pooled Frechet quantities...\n")
  s_in_plus_weights_for_pooled <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_variance_to_use, actual_kernel_function_to_pass, "+", n=n_total)
  s_in_minus_weights_for_pooled <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_variance_to_use, actual_kernel_function_to_pass, "-", n=n_total)
  s_in_pooled_weights <- 0.5 * (s_in_plus_weights_for_pooled + s_in_minus_weights_for_pooled) # Averaged weights
  active_indices_pooled <- which(s_in_pooled_weights != 0)
  
  V_hat_pooled <- NA_real_ ; l_hat_pooled_est <- NULL
  if (length(active_indices_pooled) > 0) {
    l_hat_pooled_est <- .estimate_weighted_frechet_mean(Y_obj, active_indices_pooled, 
                                                        s_in_pooled_weights[active_indices_pooled], 
                                                        metric_space_type, frechet_optns)
    if (!is.null(l_hat_pooled_est)) {
      d_sq_vals_pooled_active <- sapply(active_indices_pooled, function(i_act) {
        dist_fun_sq_calculator(Y_obj[[i_act]], l_hat_pooled_est, frechet_optns)
      })
      sum_w_pool_signed <- sum(s_in_pooled_weights[active_indices_pooled]) # Use signed sum for mean-like quantity
      if (abs(sum_w_pool_signed) < 1e-12) { # Handle case where weights nearly cancel
        warning("Sum of pooled weights is near zero. V_hat_pooled might be unstable.")
        # Fallback to unweighted average of d_sq_vals_pooled_active if sum_w_pool_signed is tiny
        V_hat_pooled <- if (length(d_sq_vals_pooled_active) > 0) mean(d_sq_vals_pooled_active, na.rm=TRUE) else NA_real_
      } else {
        V_hat_pooled <- sum(s_in_pooled_weights[active_indices_pooled] * d_sq_vals_pooled_active, na.rm = TRUE) / sum_w_pool_signed
      }
    } else {
      return(list(Tn = NA, p_value = NA, error = "Pooled Frechet mean estimation failed.",
                  h_mean_cv_selected = h_mean_cv_selected, h_variance_used = h_variance_to_use,
                  debug_info_plus = res_plus$debug_info, debug_info_minus = res_minus$debug_info))
    }
  } else {
    return(list(Tn = NA, p_value = NA, error = "No active weights for pooled estimation.",
                h_mean_cv_selected = h_mean_cv_selected, h_variance_used = h_variance_to_use,
                debug_info_plus = res_plus$debug_info, debug_info_minus = res_minus$debug_info))
  }
  if(is.na(V_hat_pooled)) V_hat_pooled <- 0 # Default if NA (e.g. if sum_w_pool_signed was zero and no active points)
  
  
  # --- Estimate f_X(c) ---
  if (verbose) cat("Estimating density f_X(c)...\n")
  h_fx_to_use <- h_fx
  if (is.null(h_fx_to_use)) {
    # Standard rule-of-thumb bandwidth for density estimation
    h_fx_temp <- tryCatch(stats::bw.nrd0(X_scalar), error = function(e) NA)
    if(is.na(h_fx_temp) || h_fx_temp <= 1e-5) { # Ensure positive and not too small
      sd_X_global <- stats::sd(X_scalar, na.rm=TRUE); if(is.na(sd_X_global) || sd_X_global < 1e-5) sd_X_global <- 1
      # Silverman's rule of thumb scaled: 1.06 * sd * n^(-1/5)
      # Using a slightly different heuristic here from original:
      h_fx_temp <- 0.9 * min(sd_X_global, stats::IQR(X_scalar, na.rm=TRUE)/1.349) * (n_total^(-1/5)) 
      if(is.na(h_fx_temp) || h_fx_temp <= 1e-5) h_fx_temp <- 0.1 # Last resort fallback
    }
    h_fx_to_use <- h_fx_temp
  }
  if (h_fx_to_use <= 1e-5) { # Final check
    sd_X_f <- stats::sd(X_scalar, na.rm=TRUE); if(is.na(sd_X_f) || sd_X_f < 1e-5) sd_X_f <- 1
    h_fx_fallback <- 0.1 * sd_X_f * (n_total^(-1/5)); if(is.na(h_fx_fallback) || h_fx_fallback <=1e-5) h_fx_fallback <- 0.01
    warning(paste("Effective h_fx was non-positive or too small. Using fallback:", signif(h_fx_fallback,3)))
    h_fx_to_use <- h_fx_fallback
  }
  
  # Ensure from/to for density are reasonable
  range_X_scalar <- range(X_scalar, na.rm = TRUE)
  density_from <- if(is.finite(range_X_scalar[1])) range_X_scalar[1] - 3 * h_fx_to_use else c_val - 5 * h_fx_to_use
  density_to   <- if(is.finite(range_X_scalar[2])) range_X_scalar[2] + 3 * h_fx_to_use else c_val + 5 * h_fx_to_use
  if (density_from >= density_to) { # Fallback if range is problematic
    density_from <- c_val - 5 * abs(c_val * 0.1 + 0.1) # Ensure some width
    density_to   <- c_val + 5 * abs(c_val * 0.1 + 0.1)
    if(density_from >= density_to) {density_from <- c_val -1; density_to <- c_val + 1;} # Last resort
  }
  
  fx_den_obj <- try(stats::density(X_scalar, bw = h_fx_to_use, kernel = kernel_fx_char, n = 1024, # n=1024 is a reasonable default
                                   from = density_from, to = density_to), silent=TRUE)
  if(inherits(fx_den_obj, "try-error")){
    warning("stats::density for f_X(c) failed. Using heuristic 1/range(X).")
    range_X_for_fx_calc <- diff(range_X_scalar)
    f_X_hat_c <- if(is.finite(range_X_for_fx_calc) && range_X_for_fx_calc > 1e-20) 1 / range_X_for_fx_calc else 1
  } else {
    f_X_hat_c <- stats::approx(fx_den_obj$x, fx_den_obj$y, xout = c_val, rule=2)$y # rule=2: use endpoint if outside range
  }
  if (is.na(f_X_hat_c) || f_X_hat_c < 1e-20) { # Floor f_X_hat_c
    warning(paste0("f_X_hat_c was NA or too small (", signif(f_X_hat_c,3), "). Setting to 1e-20."))
    f_X_hat_c <- 1e-20
  }
  
  
  # --- Kernel Constant S_K ---
  # (This section for S_K calculation appears correct as per the paper's likely intentions)
  integrate_K_moment_internal_main <- function(K_fun_char_internal_sk, k_pow, j_pow) {
    K_to_integrate_fun_sk <- get("kerFctn", envir = asNamespace("frechet"))(K_fun_char_internal_sk)
    integrand <- function(u_vec) { sapply(u_vec, function(u) { if (u < 0) return(0); K_to_integrate_fun_sk(u)^k_pow * u^j_pow }) }
    limit_val_upper <- if (K_fun_char_internal_sk %in% c("rect", "epan", "quar")) 1 else Inf # Integration limit for compact kernels
    # Try integration, with fallback for infinite upper limit
    res <- try(stats::integrate(integrand, lower = 0, upper = limit_val_upper, subdivisions = 2000, rel.tol = 1e-6)$value, silent=TRUE)
    if ((inherits(res, "try-error") || is.na(res)) && is.infinite(limit_val_upper)) { # Fallback for Gaussian-like kernels
      res <- try(stats::integrate(integrand, lower = 0, upper = 8, subdivisions = 2000, rel.tol = 1e-6)$value, silent=TRUE) # Approx integral
    }
    if (inherits(res, "try-error") || is.na(res)) NA else res
  }
  K_p_10 <- integrate_K_moment_internal_main(kernel_frechet_char, 1, 0) # mu_0(K_+)
  K_p_11 <- integrate_K_moment_internal_main(kernel_frechet_char, 1, 1) # mu_1(K_+)
  K_p_12 <- integrate_K_moment_internal_main(kernel_frechet_char, 1, 2) # mu_2(K_+)
  
  S_K_num_integrand_fun_main <- function(u_vec) {
    K_val_fun_S_sk <- get("kerFctn", envir = asNamespace("frechet"))(kernel_frechet_char)
    sapply(u_vec, function(u_s) {
      if (u_s < 0) return(0);
      K_val_at_u_sk <- K_val_fun_S_sk(u_s);
      # This is ( K_2 - u K_1 )^2 K(u)^2 from paper, where K_j = mu_j(K_+)
      return(((K_p_12 - u_s * K_p_11)^2) * (K_val_at_u_sk^2)) 
    })
  }
  limit_S_K_num_upper <- if (kernel_frechet_char %in% c("rect", "epan", "quar")) 1 else Inf
  S_num_val <- try(stats::integrate(S_K_num_integrand_fun_main, lower = 0, upper = limit_S_K_num_upper, subdivisions = 2000, rel.tol = 1e-6)$value, silent=TRUE)
  if ((inherits(S_num_val, "try-error") || is.na(S_num_val)) && is.infinite(limit_S_K_num_upper)) {
    S_num_val <- try(stats::integrate(S_K_num_integrand_fun_main, lower = 0, upper = 8, subdivisions = 2000, rel.tol = 1e-6)$value, silent=TRUE)
  }
  
  # Denominator is (K_2 K_0 - K_1^2)^2
  S_den_val_paren_sq <- if(anyNA(c(K_p_10,K_p_11,K_p_12))) NA else (K_p_12 * K_p_10 - K_p_11^2)^2
  S_K <- if(anyNA(c(S_num_val, S_den_val_paren_sq)) || is.na(S_den_val_paren_sq) || abs(S_den_val_paren_sq) < 1e-20) NA else S_num_val / S_den_val_paren_sq
  
  if(!is.na(S_K) && (!is.finite(S_K) || S_K <=0)) {
    warning(paste0("S_K calculation was non-finite or non-positive (", signif(S_K,3), "). Setting to NA."))
    S_K <- NA
  }
  if(is.na(S_K) && verbose) cat("S_K is NA. Check kernel constants: K_p_10=", K_p_10, "K_p_11=", K_p_11, "K_p_12=", K_p_12, "S_num=", S_num_val, "S_den_sq=", S_den_val_paren_sq, "\n")
  
  
  # --- Test Statistic Components ---
  F_n <- V_hat_pooled - 0.5 * (res_plus$V_hat + res_minus$V_hat)
  U_n_numerator <- (res_plus$V_hat - res_minus$V_hat)^2
  
  sigma_V_plus_for_U_den <- res_plus$sigma_V_sq_hat
  sigma_V_minus_for_U_den <- res_minus$sigma_V_sq_hat
  
  MIN_DENOM_THRESH_Tn <- 1e-25 # Threshold for denominators to avoid Inf/NaN
  
  # Denominator for U_n term. Note: S_K is for one-sided kernel.
  # For symmetric kernels, S_K^+ = S_K^-, so S_K is used for both.
  # If kernel is asymmetric, need S_K^+ and S_K^-. Here, S_K is assumed for symmetric.
  D_n_original_U <- if(is.na(S_K) || is.na(f_X_hat_c) || is.na(sigma_V_plus_for_U_den) || is.na(sigma_V_minus_for_U_den) || f_X_hat_c < MIN_DENOM_THRESH_Tn) {
    NA 
  } else {
    (S_K * sigma_V_plus_for_U_den / f_X_hat_c) + (S_K * sigma_V_minus_for_U_den / f_X_hat_c)
  }
  
  U_n <- if (is.na(D_n_original_U) || abs(D_n_original_U) < MIN_DENOM_THRESH_Tn) {
    if (!is.na(U_n_numerator) && abs(U_n_numerator) < MIN_DENOM_THRESH_Tn * 1e-3 && !is.na(D_n_original_U) && D_n_original_U != 0) 0 else NA
  } else { 
    U_n_numerator / D_n_original_U 
  }
  
  if (is.na(U_n) && !is.na(U_n_numerator) && abs(U_n_numerator) > MIN_DENOM_THRESH_Tn && 
      !is.na(D_n_original_U) && abs(D_n_original_U) < MIN_DENOM_THRESH_Tn) {
    U_n <- Inf # Numerator is substantial, denominator is zero -> Inf
  }
  
  T_n_term1_U <- if(is.na(U_n)) NA else n_total * h_variance_to_use * U_n
  
  # Denominator for F_n^2 term
  # lambda_hat definition based on sum of (absolute) weights for plus and minus sides for pooled estimation
  # sum of s_in_plus_weights_for_pooled and s_in_minus_weights_for_pooled might be more direct
  # if these are the "effective N" for each side under pooled estimation using h_variance_to_use
  
  # Using sum of absolute values of the raw one-sided weights used for pooled estimation.
  # These weights (s_in_plus/minus_weights_for_pooled) were calculated using h_variance_to_use.
  sum_abs_w_plus_for_pooled  <- sum(abs(s_in_plus_weights_for_pooled), na.rm = TRUE)
  sum_abs_w_minus_for_pooled <- sum(abs(s_in_minus_weights_for_pooled), na.rm = TRUE)
  
  lambda_hat <- if ((sum_abs_w_plus_for_pooled + sum_abs_w_minus_for_pooled) > 1e-12) {
    sum_abs_w_plus_for_pooled / (sum_abs_w_plus_for_pooled + sum_abs_w_minus_for_pooled)
  } else {
    0.5 # Default if total weight sum is tiny
  }
  
  D_n_for_Fn_sq <- if(is.na(S_K) || is.na(f_X_hat_c) || is.na(sigma_V_plus_for_U_den) || is.na(sigma_V_minus_for_U_den) || 
                      is.na(lambda_hat) || f_X_hat_c < MIN_DENOM_THRESH_Tn || lambda_hat * (1 - lambda_hat) < 1e-12 ) { # Added check for lambda_hat
    NA
  } else {
    lambda_hat * (1 - lambda_hat) * (sigma_V_plus_for_U_den + sigma_V_minus_for_U_den) * (S_K / f_X_hat_c)
  }
  
  T_n_term2_F_sq <- if (is.na(D_n_for_Fn_sq) || abs(D_n_for_Fn_sq) < MIN_DENOM_THRESH_Tn) {
    if (!is.na(F_n) && F_n^2 < MIN_DENOM_THRESH_Tn * 1e-3) { # If F_n is also tiny
      0
    } else if (!is.na(F_n) && F_n^2 >= MIN_DENOM_THRESH_Tn * 1e-3) { # If F_n is not tiny but denom is
      Inf
    } else { # F_n is NA or denom is NA
      NA
    }
  } else {
    (n_total * h_variance_to_use * F_n^2) / D_n_for_Fn_sq
  }
  
  T_n <- if(is.na(T_n_term1_U) || is.na(T_n_term2_F_sq)) NA else T_n_term1_U + T_n_term2_F_sq
  if (!is.na(T_n) && is.infinite(T_n)) T_n <- 1e15 # Cap large values
  
  p_value <- if (is.na(T_n) || T_n < 0) NA else {stats::pchisq(T_n, df = 1, lower.tail = FALSE)}
  
  # --- Prepare Return List ---
  current_error_msg <- NULL
  # Collect debug info from one-sided estimations
  if(!is.null(res_plus$debug_info) && is.character(res_plus$debug_info)) current_error_msg <- paste("Plus side:",res_plus$debug_info)
  if(!is.null(res_minus$debug_info) && is.character(res_minus$debug_info)){
    if(is.null(current_error_msg)) current_error_msg <- paste("Minus side:", res_minus$debug_info)
    else current_error_msg <- paste(current_error_msg, "; Minus side:", res_minus$debug_info)
  }
  
  # If Tn is NA, try to provide more specific reasons
  if(is.na(T_n) && is.null(current_error_msg)) {
    components_status <- c(V_plus=res_plus$V_hat, V_minus=res_minus$V_hat, V_pooled=V_hat_pooled,
                           sigV_sq_plus=sigma_V_plus_for_U_den, sigV_sq_minus=sigma_V_minus_for_U_den,
                           fx_c=f_X_hat_c, SK_const=S_K, 
                           D_U_den=D_n_original_U, D_F_den=D_n_for_Fn_sq,
                           lambda_h=lambda_hat)
    na_components <- names(components_status)[is.na(components_status)]
    zero_den_U <- if(!is.na(D_n_original_U) && abs(D_n_original_U) < MIN_DENOM_THRESH_Tn) "D_n_original_U_near_zero" else NULL
    zero_den_F <- if(!is.na(D_n_for_Fn_sq) && abs(D_n_for_Fn_sq) < MIN_DENOM_THRESH_Tn) "D_n_for_Fn_sq_near_zero" else NULL
    
    problem_flags <- c(if(length(na_components)>0) paste("NA_components:", paste(na_components, collapse=", ")), zero_den_U, zero_den_F)
    problem_flags <- problem_flags[!is.null(problem_flags)] 
    
    if(length(problem_flags) > 0) {
      current_error_msg <- paste("Tn is NA due to:", paste(problem_flags, collapse="; "))
    } else {
      current_error_msg <- "Tn calculation resulted in NA for unspecified reason (all components seemed valid before Tn calculation)."
    }
  }
  
#   # The component dump print statement (kept from original)
#   cat(sprintf(
#     "\n--- component dump (frechesTest) ---\n\
# V+=%.4e  V-=%.4e  Vp=%.4e\n\
# sigVsq+=%.4e sigVsq-=%.4e\n\
# Fn=%.4e  U_num=%.4e  U_den=%.4e  F_den=%.4e\n\
# SK=%.4e  fx_c=%.4e  h_var=%.4e  n=%d lambda_hat=%.3f\n\
# Tn1(U)=%.4e Tn2(F)=%.4e Total Tn=%.4e\n",
#     res_plus$V_hat, res_minus$V_hat, V_hat_pooled,
#     sigma_V_plus_for_U_den, sigma_V_minus_for_U_den,
#     F_n, U_n_numerator, D_n_original_U, D_n_for_Fn_sq,
#     S_K, f_X_hat_c, h_variance_to_use, n_total, lambda_hat,
#     T_n_term1_U, T_n_term2_F_sq, T_n
#   ))
  
  return(list(Tn = T_n, p_value = p_value, error = current_error_msg,
              V_hat_plus = res_plus$V_hat, V_hat_minus = res_minus$V_hat, V_hat_pooled = V_hat_pooled,
              sigma_V_sq_hat_plus = sigma_V_plus_for_U_den, sigma_V_sq_hat_minus = sigma_V_minus_for_U_den,
              f_X_hat_c = f_X_hat_c, S_K = S_K, lambda_hat = lambda_hat,
              F_n = F_n, U_n = U_n,
              l_hat_plus = res_plus$l_hat, l_hat_minus = res_minus$l_hat, l_hat_pooled = l_hat_pooled_est,
              kernel_constants = list(K_p_10=K_p_10, K_p_11=K_p_11, K_p_12=K_p_12, S_K_num = S_num_val, S_K_den_paren_sq=S_den_val_paren_sq),
              num_active_plus = res_plus$num_active_weights,
              num_active_minus = res_minus$num_active_weights,
              num_active_pooled = length(active_indices_pooled),
              h_mean_cv_selected = h_mean_cv_selected,
              h_variance_used = h_variance_to_use,
              D_n_U_denominator = D_n_original_U, 
              D_n_F_denominator = D_n_for_Fn_sq 
  ))
}