# FrechetJumpTest_helpers.R
#
# Helper functions for Frechet Jump Test.
# These functions are primarily used internally by the main Frechet Jump Test procedures.

# NULL-coalescing operator (if R version < 4.0.0)
# Provides a default value if 'a' is NULL.
`%||%` <- function(a, b) {
  if (is.null(a)) {
    b
  } else {
    a
  }
}

# Projects a target matrix B onto the space of valid Laplacian matrices using OSQP.
# Handles both directed and undirected Laplacians.
# Caches the OSQP model structure for efficiency if M_nodes and directed_status are the same.
.project_to_laplacian_space_osqp <- function(target_matrix_B,
                                             M_nodes,
                                             W_upper_edge_abs_val = 1.0,
                                             osqp_model_cache = NULL, # Expected: list(model=osqp_obj, M_nodes=.., directed_status=..)
                                             osqp_settings = NULL,
                                             is_directed = FALSE) { # Indicates if constraints for directed Laplacian should be used
  
  q_osqp_target <- -as.vector(target_matrix_B) # Target for OSQP: -B in min ||L - B||^2_F
  model_to_solve <- NULL
  
  # --- Cache Handling ---
  # Determine if the existing cache is valid for the current request (M_nodes and directed status).
  # The frechesTest function is responsible for creating/passing the cache object
  # that includes the 'directed_status' it was built for.
  
  # This helper assumes osqp_model_cache IS the model object if valid,
  # otherwise it builds a new one. The calling function (frechesTest)
  # handles storing this newly built model back into frechet_optns if needed.
  
  rebuild_this_call <- TRUE # Default to rebuild if cache is not perfectly matching
  if (!is.null(osqp_model_cache) &&
      !is.null(osqp_model_cache$model) &&
      M_nodes == osqp_model_cache$M_nodes &&
      !is.null(osqp_model_cache$directed_status) &&
      is_directed == osqp_model_cache$directed_status) {
    rebuild_this_call <- FALSE
    model_to_solve <- osqp_model_cache$model
  }
  
  # --- OSQP Model Building (if necessary) ---
  if (rebuild_this_call) {
    # Define the P matrix for the quadratic program (||L - B||_F^2 means P = I)
    P_osqp_sparse <- Matrix::Matrix(diag(M_nodes^2), sparse = TRUE)
    
    # Determine number of constraints based on directed status
    num_off_diag_total_possible <- M_nodes * (M_nodes - 1)
    num_off_diag_unique_pairs <- num_off_diag_total_possible / 2
    
    num_constraints_sym <- if (!is_directed && M_nodes > 1) num_off_diag_unique_pairs else 0
    num_constraints_rowsum <- M_nodes
    num_constraints_offdiag_bounds <- if (M_nodes > 1) {
      if (is_directed) num_off_diag_total_possible else num_off_diag_unique_pairs
    } else {
      0
    }
    num_total_constraints <- num_constraints_sym + num_constraints_rowsum + num_constraints_offdiag_bounds
    
    # Initialize constraint structures for sparse matrix A
    # Using lists to build sparse matrix components efficiently before creating the sparseMatrix object.
    row_indices   <- integer(0)
    col_indices   <- integer(0)
    values        <- numeric(0)
    l_bounds_list <- numeric(0)
    u_bounds_list <- numeric(0)
    current_max_row <- 0
    
    # Constraint 1: Symmetry (L_rc - L_cr = 0) --- ONLY IF UNDIRECTED
    if (!is_directed && num_constraints_sym > 0) {
      for (c_node in 2:M_nodes) {
        for (r_node in 1:(c_node - 1)) { # Iterate over upper triangle (r_node < c_node)
          current_max_row <- current_max_row + 1
          idx_rc <- (c_node - 1) * M_nodes + r_node # Element L_r_node,c_node (column-major index)
          idx_cr <- (r_node - 1) * M_nodes + c_node # Element L_c_node,r_node
          
          row_indices   <- c(row_indices, current_max_row, current_max_row)
          col_indices   <- c(col_indices, idx_rc, idx_cr)
          values        <- c(values, 1, -1) # L_rc - L_cr
          
          l_bounds_list <- c(l_bounds_list, 0)  # = 0
          u_bounds_list <- c(u_bounds_list, 0)  # = 0
        }
      }
    }
    
    # Constraint 2: Zero Row Sums (sum_j L_ij = 0 for all i)
    for (r_node in 1:M_nodes) {
      current_max_row <- current_max_row + 1
      # Indices for all elements in row r_node: L_r_node,1, L_r_node,2, ..., L_r_node,M_nodes
      indices_in_row <- ((1:M_nodes) - 1) * M_nodes + r_node 
      
      row_indices   <- c(row_indices, rep(current_max_row, M_nodes))
      col_indices   <- c(col_indices, indices_in_row)
      values        <- c(values, rep(1, M_nodes)) # Sum of elements in row r_node
      
      l_bounds_list <- c(l_bounds_list, 0)  # = 0
      u_bounds_list <- c(u_bounds_list, 0)  # = 0
    }
    
    # Constraint 3: Off-diagonal Bounds (-W <= L_ij <= 0 for i != j)
    if (num_constraints_offdiag_bounds > 0) {
      if (!is_directed) { # For UNDIRECTED, only need to constrain unique off-diagonal pairs (e.g., upper triangle)
        for (c_node in 2:M_nodes) {
          for (r_node in 1:(c_node - 1)) { # Iterate over upper triangle (r_node < c_node)
            current_max_row <- current_max_row + 1
            idx_rc <- (c_node - 1) * M_nodes + r_node # Element L_r_node,c_node
            
            row_indices   <- c(row_indices, current_max_row)
            col_indices   <- c(col_indices, idx_rc)
            values        <- c(values, 1) # L_r_node,c_node itself
            
            l_bounds_list <- c(l_bounds_list, -W_upper_edge_abs_val) # -W <= L_rc
            u_bounds_list <- c(u_bounds_list, 0)                     # L_rc <= 0
          }
        }
      } else { # For DIRECTED, constrain ALL M_nodes * (M_nodes - 1) off-diagonal elements
        for (r_constr_node in 1:M_nodes) {
          for (c_constr_node in 1:M_nodes) {
            if (r_constr_node == c_constr_node) next # Skip diagonal elements
            
            current_max_row <- current_max_row + 1
            idx_rc <- (c_constr_node - 1) * M_nodes + r_constr_node # Element L_r_constr_node,c_constr_node
            
            row_indices   <- c(row_indices, current_max_row)
            col_indices   <- c(col_indices, idx_rc)
            values        <- c(values, 1) # L_r_constr_node,c_constr_node itself
            
            l_bounds_list <- c(l_bounds_list, -W_upper_edge_abs_val) # -W <= L_rc
            u_bounds_list <- c(u_bounds_list, 0)                     # L_rc <= 0
          }
        }
      }
    }
    
    # Assemble the constraint matrix A and bounds vectors l, u
    if (num_total_constraints > 0) {
      const_matrix_A_sparse <- Matrix::sparseMatrix(i = row_indices, 
                                                    j = col_indices, 
                                                    x = values,
                                                    dims = c(current_max_row, M_nodes^2))
      l_bounds_vec <- l_bounds_list
      u_bounds_vec <- u_bounds_list
    } else { 
      # This case should be rare, e.g. M_nodes = 1 if row sum wasn't added, or if no constraints generated.
      # For M_nodes = 1, only L11=0 due to row sum, so num_total_constraints is typically 1.
      const_matrix_A_sparse <- Matrix::Matrix(0, nrow = 0, ncol = M_nodes^2, sparse = TRUE)
      l_bounds_vec <- numeric(0)
      u_bounds_vec <- numeric(0)
    }
    
    # Special handling for M_nodes = 1 if no constraints were added (e.g. row sum logic error for M_nodes=1).
    # This ensures L11=0 constraint is present.
    if (M_nodes == 1 && num_total_constraints == 0) { 
      const_matrix_A_sparse <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 1))
      l_bounds_vec <- 0
      u_bounds_vec <- 0
    }
    
    # Setup and initialize the OSQP problem
    current_osqp_settings <- osqp_settings %||% osqp::osqpSettings(
      verbose = FALSE, 
      eps_abs = 1e-8, 
      eps_rel = 1e-8, 
      max_iter = 20000, 
      polish = TRUE
    )
    
    model_to_solve <- osqp::osqp(
      P = P_osqp_sparse, 
      q = q_osqp_target, # q_osqp_target is based on current data matrix B
      A = const_matrix_A_sparse, 
      l = l_bounds_vec, 
      u = u_bounds_vec,
      pars = current_osqp_settings
    )
    # Note: The calling function (e.g., frechesTest) should update its cache with:
    # list(model = model_to_solve, M_nodes = M_nodes, directed_status = is_directed, cache_id = new_cache_id)
    
  } else { # Cache is valid, model_to_solve is already set from cache
    # Update only the q vector (target) for the existing model
    model_to_solve$Update(q = q_osqp_target)
  }
  
  # --- Solve the OSQP Problem ---
  solution <- model_to_solve$Solve()
  
  # --- Fallback Mechanism (if OSQP fails) ---
  # Check solution status and apply fallback constraints if not optimal
  if (solution$info$status_val != 1 && solution$info$status_val != 2 ) { # 1: solved, 2: solved inaccurate
    warning(paste("OSQP projection for", if (is_directed) "directed" else "undirected",
                  "Laplacian did not solve optimally. Status:", solution$info$status,
                  "(Val:", solution$info$status_val, "). Applying fallback constraints."))
    
    projected_L_matrix <- target_matrix_B # Start with the original target matrix
    
    # Symmetrize if an undirected graph is expected
    if (!is_directed) { 
      projected_L_matrix <- (target_matrix_B + t(target_matrix_B)) / 2
    }
    
    # Apply off-diagonal constraints: L_ij <= 0 and L_ij >= -W
    for (r_idx in 1:M_nodes) { 
      for (c_idx in 1:M_nodes) {
        if (r_idx != c_idx) {
          projected_L_matrix[r_idx, c_idx] <- min(0, max(-W_upper_edge_abs_val, projected_L_matrix[r_idx, c_idx]))
        }
      }
    }
    
    # Enforce zero row sums by adjusting diagonal elements
    diag(projected_L_matrix) <- 0 # Temporarily zero out diagonal to sum off-diagonals
    diag(projected_L_matrix) <- -rowSums(projected_L_matrix) # Set L_ii = -sum_{j!=i} L_ij
    
    # For undirected graphs, re-symmetrize if off-diagonal constraints might have broken symmetry,
    # and re-ensure diagonal correctly reflects sum of new off-diagonals.
    if (!is_directed) {
      projected_L_matrix <- (projected_L_matrix + t(projected_L_matrix)) / 2
      diag(projected_L_matrix) <- 0 # Zero out again before recalculating row sums
      diag(projected_L_matrix) <- -rowSums(projected_L_matrix)
    }
    return(projected_L_matrix)
  }
  
  # --- Successful OSQP Solution ---
  # If solution is optimal, reshape the solution vector to a matrix
  projected_L_vector <- solution$x
  projected_L_matrix <- matrix(projected_L_vector, nrow = M_nodes, ncol = M_nodes, byrow = FALSE) # R matrices are column-major
  
  return(projected_L_matrix)
}


# Converts various observation types (raw data, pre-computed quantiles, histograms)
# into quantile functions on a target support grid (qSup_target).
.get_quantiles_for_obs_jump_test <- function(obs, qSup_target, den_opts) {
  # This function might be called less for its expensive branches if pre-computation is used effectively.
  
  # Case 1: 'obs' is already a pre-computed quantile object
  if (is.list(obs) && !is.null(obs$type) && obs$type == "quantile" && 
      !is.null(obs$q) && !is.null(obs$qSup)) {
    
    # If quantiles are already on the target qSup grid
    if (length(obs$qSup) == length(qSup_target) && all(abs(obs$qSup - qSup_target) < 1e-9)) {
      return(obs$q) 
    } else {
      # Interpolate if on a different qSup grid
      # (Ideally avoided by ensuring consistent qSup throughout the main process for performance)
      warning("Interpolating pre-computed quantiles to a new qSup_target. Ensure qSup is consistent for best performance.")
      return(stats::approx(x = obs$qSup, y = obs$q, xout = qSup_target, rule = 2)$y)
    }
    
    # Case 2: 'obs' is raw numeric data
  } else if (is.numeric(obs)) {
    # Perform density estimation from raw data
    current_den_opts <- den_opts %||% list()
    
    # Rename options for frechet::CreateDensity if old names are used (for backward compatibility)
    if (!is.null(current_den_opts$kernelDen)) { names(current_den_opts)[names(current_den_opts) == "kernelDen"] <- "kernel" }
    if (!is.null(current_den_opts$bwDen))     { names(current_den_opts)[names(current_den_opts) == "bwDen"]     <- "userBwMu" }
    if (!is.null(current_den_opts$ndSup))     { names(current_den_opts)[names(current_den_opts) == "ndSup"]     <- "nRegGrid" }
    if (!is.null(current_den_opts$dSup))      { names(current_den_opts)[names(current_den_opts) == "dSup"]      <- "outputGrid" }
    
    density_obj <- frechet::CreateDensity(y = obs, optns = current_den_opts)
    quantiles   <- fdadensity::dens2quantile(dens = density_obj$y, dSup = density_obj$x, qSup = qSup_target)
    return(quantiles)
    
    # Case 3: 'obs' is a histogram object (list with mids, counts, breaks)
  } else if (is.list(obs) && !is.null(obs$mids) && !is.null(obs$counts) && !is.null(obs$breaks)) {
    # Perform density estimation from histogram data
    current_den_opts <- den_opts %||% list()
    
    # Rename options for frechet::CreateDensity if old names are used
    if (!is.null(current_den_opts$kernelDen)) { names(current_den_opts)[names(current_den_opts) == "kernelDen"] <- "kernel" }
    if (!is.null(current_den_opts$bwDen))     { names(current_den_opts)[names(current_den_opts) == "bwDen"]     <- "userBwMu" }
    if (!is.null(current_den_opts$ndSup))     { names(current_den_opts)[names(current_den_opts) == "ndSup"]     <- "nRegGrid" }
    if (!is.null(current_den_opts$dSup))      { names(current_den_opts)[names(current_den_opts) == "dSup"]      <- "outputGrid" }
    
    density_obj <- frechet::CreateDensity(histogram = obs, optns = current_den_opts)
    quantiles   <- fdadensity::dens2quantile(dens = density_obj$y, dSup = density_obj$x, qSup = qSup_target)
    return(quantiles)
    
    # Case 4: Unsupported observation type
  } else {
    stop("Unsupported observation type for .get_quantiles_for_obs_jump_test. Expected numeric vector, specific quantile list, or histogram list.")
  }
}

# Calculates one-sided local polynomial regression weights.
# These weights are used to estimate one-sided Frechet means/variances around an evaluation point 'c_val'.
# 'n' is the number of observations on the specific side being considered.
.calculate_one_sided_locpoly_weights <- function(x_obs, c_val, h_val, K_fun, side, n) {
  # Calculate kernel values scaled by bandwidth
  Kh_X_minus_c_all <- K_fun((x_obs - c_val) / h_val) / h_val
  
  mu_hat_j <- numeric(3) # To store moments mu_0, mu_1, mu_2 for local polynomial regression
  
  # Determine which observations are on the specified 'side' ("+" or "-") of c_val
  indicator <- if (side == "+") {
    x_obs >= c_val
  } else { # side == "-"
    x_obs < c_val
  }
  active_indices_side <- which(indicator)
  
  # Check if enough points are available on the side for stable estimation
  if (length(active_indices_side) < 2) {
    warning(paste("Too few points on side", side, 
                  "for local polynomial regression (need at least 2). Weights will be zero."))
    return(rep(0, length(x_obs))) # Return zero weights for all original observations
  }
  
  # Use only data on the active side for moment calculations
  x_obs_side          <- x_obs[active_indices_side]
  Kh_X_minus_c_side   <- Kh_X_minus_c_all[active_indices_side]
  
  # Calculate moments (S_j terms in Fan & Gijbels notation, or mu_hat here)
  for (j_idx in 0:2) { # For local linear, we need up to j=2 for the weight formula
    mu_hat_j[j_idx + 1] <- sum(Kh_X_minus_c_side * (x_obs_side - c_val)^j_idx)
  }
  
  # Calculate sigma_0_sq_hat (denominator term for weights, related to variance of estimator)
  # This is (S_0 * S_2 - S_1^2) in some notations
  sigma_0_sq_hat <- mu_hat_j[1] * mu_hat_j[3] - mu_hat_j[2]^2
  
  s_in <- rep(0, length(x_obs)) # Initialize weights to zero for all original observations
  
  # Check for near-zero denominator to avoid instability
  if (abs(sigma_0_sq_hat) < 1e-20) { 
    warning(paste("sigma_0_sq_hat is close to zero for side", side, 
                  "(value:", sigma_0_sq_hat, "). Weights might be unstable or zero."))
    # s_in remains zero in this case
  } else {
    # Calculate weights for the observations on the active side
    # Formula for weights for estimating f(c) (the intercept term alpha_0)
    s_in[active_indices_side] <- (Kh_X_minus_c_side / sigma_0_sq_hat) * 
      (mu_hat_j[3] - mu_hat_j[2] * (x_obs_side - c_val))
  }
  
  return(s_in)
}

# Estimates the weighted Frechet mean for various metric spaces.
# Y_obj_list: list of objects in the metric space.
# active_idx: indices of Y_obj_list to consider for the mean.
# s_in_weights_active: weights corresponding to the active_idx.
.estimate_weighted_frechet_mean <- function(Y_obj_list, 
                                            active_idx, 
                                            s_in_weights_active,
                                            metric_space_type, 
                                            frechet_optns) {
  
  if (length(active_idx) == 0) {
    warning("No active samples for weighted Frechet mean estimation.")
    return(NULL)
  }
  
  Y_active        <- Y_obj_list[active_idx] # Subset of objects for mean calculation
  l_hat_mean      <- NULL                   # Initialize result
  sum_s_in_active <- sum(s_in_weights_active)
  
  # Determine if fallback to unweighted mean is needed (e.g., if sum of weights is effectively zero)
  # Sphere space has its own robust handling for the initial guess in trust region method.
  use_fallback_unweighted_mean <- abs(sum_s_in_active) < 1e-9 && metric_space_type != "sphere"
  
  # --- Fallback Logic (if sum of weights is near zero) ---
  if (use_fallback_unweighted_mean) {
    warning("Sum of active weights is close to zero. Using unweighted Frechet mean of active samples as fallback.")
    
    if (metric_space_type == "density") {
      qSup           <- frechet_optns$qSup %||% stop("qSup missing for density fallback mean")
      den_opts_cfd   <- frechet_optns$den_opts_for_create_density %||% list()
      
      # Convert active Y objects to quantiles if they are not already
      if (is.list(Y_active[[1]]) && !is.null(Y_active[[1]]$type) && Y_active[[1]]$type == "quantile") {
        qin_active_fallback <- do.call(rbind, lapply(Y_active, `[[`, "q"))
      } else {
        qin_active_fallback <- t(sapply(Y_active, .get_quantiles_for_obs_jump_test, 
                                        qSup_target = qSup, den_opts = den_opts_cfd))
      }
      l_hat_q_unprojected <- colMeans(qin_active_fallback, na.rm = TRUE) # Unweighted mean of quantiles
      
    } else if (metric_space_type %in% c("covariance", "correlation", "network")) {
      l_hat_mat_unproj <- Reduce("+", Y_active) / length(Y_active) # Unweighted mean of matrices
      
    } else {
      # Sphere type is handled by trust method which can start from an arbitrary point if y0 is poor.
      stop("Fallback for zero sum of weights not implemented for this space type. Sphere is handled by trust.")
    }
    
    # --- Standard Weighted Frechet Mean Calculation ---
  } else { 
    
    # --- DENSITY SPACE ---
    if (metric_space_type == "density") {
      qSup         <- frechet_optns$qSup %||% stop("qSup missing for density mean")
      den_opts_cfd <- frechet_optns$den_opts_for_create_density %||% list()
      
      # Convert all active observations to quantiles on the target qSup grid
      qin_matrix_list <- lapply(Y_active, function(obs_item) {
        .get_quantiles_for_obs_jump_test(obs_item, qSup_target = qSup, den_opts = den_opts_cfd)
      })
      qin_active_matrix   <- do.call(rbind, qin_matrix_list)
      # Weighted average of quantile functions (unprojected)
      l_hat_q_unprojected <- as.vector(crossprod(qin_active_matrix, s_in_weights_active) / sum_s_in_active)
      
      # --- COVARIANCE, CORRELATION, NETWORK SPACES (Matrices) ---
    } else if (metric_space_type %in% c("covariance", "correlation", "network")) {
      M_active_list <- Y_active
      metric        <- frechet_optns$metric %||% "frobenius"
      alpha         <- frechet_optns$alpha %||% 1 # Power for power Euclidean metric (if applicable)
      
      # Calculate unprojected mean based on the specified metric for matrices
      if (metric %in% c("frobenius", "power")) {
        # Case 1: Frobenius for network, or Frobenius/Power(alpha=1) for cov/corr (standard Euclidean mean)
        if (metric_space_type == "network" || 
            (metric_space_type %in% c("covariance", "correlation") && metric == "frobenius") ||
            (metric_space_type %in% c("covariance", "correlation") && alpha == 1)) {
          
          l_hat_mat_unproj <- Reduce("+", mapply("*", M_active_list, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
          
          # Case 2: Power Euclidean metric with alpha != 1 for cov/corr
        } else { 
          if (alpha == 0) { # Log-Euclidean metric
            log_M_active <- lapply(M_active_list, function(mat) {
              eig          <- eigen(mat, symmetric = TRUE)
              P            <- eig$vectors
              Lambda_log   <- diag(log(pmax(1e-30, eig$values)), nrow = length(eig$values)) # pmax for numerical stability
              P %*% Lambda_log %*% t(P) 
            })
            mean_log_M       <- Reduce("+", mapply("*", log_M_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
            eig_mean_log     <- eigen(mean_log_M, symmetric = TRUE)
            l_hat_mat_unproj <- eig_mean_log$vectors %*% 
              diag(exp(eig_mean_log$values), nrow = length(eig_mean_log$values)) %*% 
              t(eig_mean_log$vectors)
          } else { # Power Euclidean with alpha != 0, 1
            M_alpha_active <- lapply(M_active_list, function(mat) {
              eig          <- eigen(mat, symmetric = TRUE)
              P            <- eig$vectors
              Lambda_alpha <- diag(pmax(0, eig$values)^alpha, nrow = length(eig$values)) # pmax for non-negative eigenvalues
              P %*% Lambda_alpha %*% t(P)
            })
            mean_M_alpha     <- Reduce("+", mapply("*", M_alpha_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
            eig_mean_alpha   <- eigen(mean_M_alpha, symmetric = TRUE)
            l_hat_mat_unproj <- eig_mean_alpha$vectors %*% 
              diag(pmax(0, eig_mean_alpha$values)^(1/alpha), nrow = length(eig_mean_alpha$values)) %*% 
              t(eig_mean_alpha$vectors)
          }
        }
      } else if (metric == "log_cholesky") { # Log-Cholesky metric
        chol_parts_active <- lapply(M_active_list, function(mat) {
          LL          <- chol(mat)
          L_part      <- LL - diag(diag(LL)) # Lower triangular part (off-diagonal)
          D_diag_part <- diag(LL)            # Diagonal elements
          list(L = L_part, D_log_diag = log(pmax(1e-30, D_diag_part))) # Log of diagonal
        })
        mean_L           <- Reduce("+", mapply("*", lapply(chol_parts_active, `[[`, "L"), s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        mean_D_log_diag  <- Reduce("+", mapply("*", lapply(chol_parts_active, `[[`, "D_log_diag"), s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        SS               <- mean_L + diag(exp(mean_D_log_diag)) # Reconstruct Cholesky factor
        l_hat_mat_unproj <- t(SS) %*% SS                       # Reconstruct matrix
        
      } else if (metric == "cholesky") { # Cholesky metric (Euclidean on Cholesky factors)
        L_chol_active    <- lapply(M_active_list, chol)
        mean_L_chol      <- Reduce("+", mapply("*", L_chol_active, s_in_weights_active, SIMPLIFY = FALSE)) / sum_s_in_active
        l_hat_mat_unproj <- t(mean_L_chol) %*% (mean_L_chol)
        
      } else { 
        stop(paste("Unsupported metric for matrix spaces:", metric)) 
      }
      
      # --- SPHERE SPACE ---
    } else if (metric_space_type == "sphere") {
      yin_active_mat <- do.call(rbind, Y_active) # Combine list of vectors into a matrix
      
      # Initial guess (y0) for Frechet mean on sphere: weighted average, then normalize
      y0_unnorm      <- colSums(yin_active_mat * s_in_weights_active, na.rm = TRUE)
      norm_y0_unnorm <- frechet:::l2norm(y0_unnorm) # Using frechet's internal L2 norm
      
      if (any(is.na(y0_unnorm)) || norm_y0_unnorm < 1e-9) {
        # Fallback initial guess if weighted sum is problematic (e.g., all weights cancel, NaNs, or zero norm)
        y0 <- yin_active_mat[1,]; # Use the first active point
        if (is.null(y0)) { return(NULL) } # Should not happen if active_idx > 0
      } else {
        y0 <- y0_unnorm / norm_y0_unnorm # Normalize to be on the sphere
      }
      
      # Objective function for trust region optimization on the sphere
      # Minimizes sum of weighted squared geodesic distances
      objFctn_sphere = function(y_sphere_candidate) {
        # Ensure candidate point is on the sphere (for safety, though trust usually respects constraints)
        if (!isTRUE(all.equal(frechet:::l2norm(y_sphere_candidate), 1))) {
          return(list(value = Inf, 
                      gradient = rep(Inf, length(y_sphere_candidate)), 
                      hessian = diag(Inf, length(y_sphere_candidate))))
        }
        
        # Calculate sum of weighted squared geodesic distances
        dists_sq_vals <- sapply(1:nrow(yin_active_mat), function(i) {
          frechet::SpheGeoDist(yin_active_mat[i,], y_sphere_candidate)^2
        })
        # f_val <- sum(s_in_weights_active * dists_sq_vals, na.rm = TRUE) # Value only needed for debug
        
        # Calculate gradient and Hessian of the objective function
        g_val_terms <- matrix(0, nrow = nrow(yin_active_mat), ncol = ncol(yin_active_mat))
        h_val_sum   <- matrix(0, ncol = ncol(yin_active_mat), nrow = ncol(yin_active_mat))
        
        for (i_s in 1:nrow(yin_active_mat)) {
          dist_val <- frechet::SpheGeoDist(yin_active_mat[i_s,], y_sphere_candidate)
          # Avoid issues at 0 or pi distances where gradient/Hessian might be ill-defined/problematic
          if (abs(dist_val) > 1e-8 && abs(dist_val - pi) > 1e-8) {
            grad_i_val_unscaled <- frechet::SpheGeoGrad(yin_active_mat[i_s,], y_sphere_candidate)
            hess_i_val_unscaled <- frechet::SpheGeoHess(yin_active_mat[i_s,], y_sphere_candidate)
            
            g_val_terms[i_s,] <- s_in_weights_active[i_s] * 2 * dist_val * grad_i_val_unscaled
            h_val_sum         <- h_val_sum + s_in_weights_active[i_s] * 2 * 
              (grad_i_val_unscaled %*% t(grad_i_val_unscaled) + dist_val * hess_i_val_unscaled)
          }
        }
        return(list(value = sum(s_in_weights_active * dists_sq_vals, na.rm = TRUE), 
                    gradient = colSums(g_val_terms, na.rm = TRUE), 
                    hessian = h_val_sum))
      }
      
      # Perform optimization using trust region method
      res_trust <- try(trust::trust(objFctn_sphere, parinit = y0, rinit = 0.1, rmax = 100, iterlim = 100), silent = TRUE)
      
      if (inherits(res_trust, "try-error") || !res_trust$converged) {
        l_hat_mean <- y0 # Fallback to initial guess if optimization fails
      } else {
        l_hat_mean <- res_trust$argument / frechet:::l2norm(res_trust$argument) # Ensure result is on sphere
      }
      
      ## --- NEW: optional reflection into the positive patch -------------
      if (isTRUE(frechet_optns$enforce_positive)) {
        if (any(l_hat_mean < 0)) {
          l_hat_mean <- abs(l_hat_mean)
          l_hat_mean <- l_hat_mean / sqrt(sum(l_hat_mean^2))
        }
      }
      
    } else { # Unsupported metric space type
      stop(paste("Unsupported metric space type:", metric_space_type))
    }
  } # End of standard/fallback weighted Frechet mean calculation
  
  # --- Final Projection/Formatting Step (for certain metric spaces) ---
  
  # For DENSITY space: Project estimated quantile function to be monotone and within bounds (if specified)
  if (metric_space_type == "density") {
    m_dim <- length(l_hat_q_unprojected) # Dimension of the quantile function vector
    
    # OSQP model caching for density projection (based on dimension and bounds)
    cache_key_density <- paste("density_proj", m_dim, 
                               if (!is.null(frechet_optns$lower)) digest::digest(frechet_optns$lower) else "NULL_L", 
                               if (!is.null(frechet_optns$upper)) digest::digest(frechet_optns$upper) else "NULL_U")
    
    model_obj_density_proj <- NULL
    if (is.list(frechet_optns) && 
        !is.null(frechet_optns$osqp_model_for_density_proj) && 
        !is.null(frechet_optns$osqp_model_for_density_proj$key) && 
        frechet_optns$osqp_model_for_density_proj$key == cache_key_density) {
      model_obj_density_proj <- frechet_optns$osqp_model_for_density_proj$model
    }
    
    if (is.null(model_obj_density_proj)) { # Build OSQP model if not cached or cache is invalid
      Pmat_qp      <- as(diag(m_dim), "sparseMatrix") # P = I for ||q_proj - q_unproj||^2
      # Monotonicity constraints: q_{k+1} - q_k >= 0  => -q_k + q_{k+1} >= 0
      A_mono_i     <- rep(1:(m_dim - 1), each = 2)
      A_mono_j     <- unlist(lapply(1:(m_dim - 1), function(k) c(k, k + 1)))
      A_mono_vals  <- rep(c(-1, 1), m_dim - 1) # Coefficients for -q_k + q_{k+1}
      Amat_qp_mono <- Matrix::sparseMatrix(i = A_mono_i, j = A_mono_j, x = A_mono_vals, dims = c(m_dim - 1, m_dim))
      
      l_bounds_mono     <- rep(0, m_dim - 1)     # Lower bound for monotonicity (>=0)
      u_bounds_mono     <- rep(Inf, m_dim - 1)   # Upper bound for monotonicity (no explicit upper)
      Final_Amat_qp     <- Amat_qp_mono
      Final_l_bounds_qp <- l_bounds_mono
      Final_u_bounds_qp <- u_bounds_mono
      
      # Optional lower bound for q_1 (first element of quantile function)
      if (!is.null(frechet_optns$lower)) {
        Final_Amat_qp     <- rbind(Final_Amat_qp, Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, m_dim)))
        Final_l_bounds_qp <- c(Final_l_bounds_qp, frechet_optns$lower)
        Final_u_bounds_qp <- c(Final_u_bounds_qp, Inf) # q_1 >= lower
      }
      # Optional upper bound for q_m (last element of quantile function)
      if (!is.null(frechet_optns$upper)) {
        Final_Amat_qp     <- rbind(Final_Amat_qp, Matrix::sparseMatrix(i = 1, j = m_dim, x = 1, dims = c(1, m_dim)))
        Final_l_bounds_qp <- c(Final_l_bounds_qp, -Inf) # q_m <= upper
        Final_u_bounds_qp <- c(Final_u_bounds_qp, frechet_optns$upper)
      }
      
      current_osqp_settings_density <- frechet_optns$osqp_settings_density_proj %||% 
        frechet_optns$osqp_settings %||% 
        osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8)
      q_placeholder_density <- rep(0, m_dim) # Placeholder for q vector in OSQP model definition
      
      model_obj_density_proj <- osqp::osqp(P = Pmat_qp, q = q_placeholder_density, 
                                           A = Final_Amat_qp, l = Final_l_bounds_qp, u = Final_u_bounds_qp, 
                                           pars = current_osqp_settings_density)
      # Store the built model in frechet_optns for caching in subsequent calls
      if (is.list(frechet_optns)) { 
        # This modification of frechet_optns might need careful handling if R's pass-by-value for lists is an issue.
        # However, if frechet_optns is an environment or list within an environment, changes persist.
        # Assuming the main function handles this correctly.
        frechet_optns[['osqp_model_for_density_proj']] <- list(model = model_obj_density_proj, key = cache_key_density)
      }
    }
    
    # Update q vector with -l_hat_q_unprojected and solve the projection QP
    model_obj_density_proj$Update(q = -l_hat_q_unprojected) 
    osqp_sol <- model_obj_density_proj$Solve()
    
    l_hat_q_projected <- if (osqp_sol$info$status_val != 1 && osqp_sol$info$status_val != 2 ) { 
      warning("OSQP for density projection failed. Returning unprojected quantiles.")
      l_hat_q_unprojected # Fallback to unprojected if OSQP fails
    } else { 
      osqp_sol$x # Projected quantiles
    }
    l_hat_mean <- list(q = l_hat_q_projected, qSup = frechet_optns$qSup, type = "quantile")
    
    # For COVARIANCE, CORRELATION, NETWORK spaces: Project to appropriate matrix space
  } else if (metric_space_type %in% c("covariance", "correlation", "network")) {
    is_network_directed_local <- frechet_optns$network_directed %||% FALSE # Default to undirected if not specified
    
    # Symmetrize for undirected networks, covariance, or correlation.
    # For directed networks, use the (potentially) asymmetric unprojected mean.
    target_matrix_for_proj <- if (metric_space_type == "network" && is_network_directed_local) {
      l_hat_mat_unproj # Use as-is for directed networks
    } else {
      (l_hat_mat_unproj + t(l_hat_mat_unproj)) / 2 # Ensure symmetry
    }
    
    if (metric_space_type == "network") {
      # Project to Laplacian space using the dedicated OSQP helper
      num_nodes_in_data <- ncol(Y_obj_list[[1]]) # Assumes all matrices in list have same dimensions
      W_for_projection  <- frechet_optns$W_laplacian_bound %||% 1.0
      
      l_hat_mean <- .project_to_laplacian_space_osqp(
        target_matrix_B      = target_matrix_for_proj,
        M_nodes              = num_nodes_in_data,
        W_upper_edge_abs_val = W_for_projection,
        osqp_model_cache     = frechet_optns$osqp_model_for_laplacian_proj, # Cache includes directed_status
        osqp_settings        = frechet_optns$osqp_settings_laplacian_proj,
        is_directed          = is_network_directed_local # Pass the directed flag
      )
      
    } else { # Covariance / Correlation: Project to Positive Semi-Definite (PSD) space
      is_corr <- (metric_space_type == "correlation")
      
      # Use Matrix::nearPD to find the nearest PSD matrix
      # target_matrix_for_proj is already symmetrized for cov/corr at this point
      l_hat_mat_pd <- try(Matrix::nearPD(target_matrix_for_proj, 
                                         corr = is_corr,           # Ensure unit diagonal if correlation
                                         ensureSymmetry = TRUE,    # Input should be symmetric, but enforce
                                         keepDiag = is_corr,       # For correlation, try to keep diagonal 1
                                         maxit = 100)$mat, 
                          silent = TRUE)
      
      if (inherits(l_hat_mat_pd, "try-error")) {
        warning(paste("nearPD failed for", metric_space_type, "mean. Using input to nearPD as fallback."))
        l_hat_mat <- target_matrix_for_proj # Fallback to the (symmetrized) unprojected mean
      } else {
        l_hat_mat <- as.matrix(l_hat_mat_pd)
      }
      
      if (is_corr) {
        diag(l_hat_mat) <- 1 # Strictly enforce unit diagonal for correlation matrices
        # Re-symmetrize after fixing diagonal, just in case nearPD didn't perfectly preserve it or it was altered.
        l_hat_mat <- (l_hat_mat + t(l_hat_mat)) / 2
      }
      l_hat_mean <- l_hat_mat
    }
  }
  # For sphere space, l_hat_mean is already set by the trust region optimization and is on the sphere.
  
  return(l_hat_mean)
}

# Estimates one-sided Frechet mean (l_hat), one-sided Frechet variance (V_hat), 
# and the variance of V_hat (sigma_V_sq_hat).
.estimate_one_sided_frechet_quantities <- function(Y_obj_list, 
                                                   X_scalar, 
                                                   metric_space_type,
                                                   c_val, 
                                                   h_frechet, 
                                                   K_frechet_fun, 
                                                   side,
                                                   frechet_optns, 
                                                   dist_fun_sq_calculator) {
  n_total <- length(Y_obj_list)
  
  # 1. Calculate one-sided local polynomial weights for observations on the specified 'side'
  s_in_weights <- .calculate_one_sided_locpoly_weights(X_scalar, c_val, h_frechet, 
                                                       K_frechet_fun, side, n_total)
  
  active_indices <- which(s_in_weights != 0) # Indices of observations with non-zero weights
  if (length(active_indices) == 0) {
    warning(paste("No data points have non-zero weight for side", side, 
                  ". Cannot estimate Frechet mean/variance."))
    return(list(l_hat = NULL, V_hat = NA, sigma_V_sq_hat = NA, 
                num_active_weights = 0, debug_info = "No active weights"))
  }
  s_in_weights_active <- s_in_weights[active_indices]
  
  # 2. Estimate the weighted Frechet mean (l_hat) using active samples and their weights
  # This call uses the .estimate_weighted_frechet_mean function, which handles various spaces and projections.
  l_hat <- .estimate_weighted_frechet_mean(Y_obj_list          = Y_obj_list, 
                                           active_idx          = active_indices, 
                                           s_in_weights_active = s_in_weights_active,
                                           metric_space_type   = metric_space_type, 
                                           frechet_optns       = frechet_optns) 
  
  if (is.null(l_hat)) {
    # If Frechet mean estimation failed (e.g., no data, or internal error in .estimate_weighted_frechet_mean)
    return(list(l_hat = NULL, V_hat = NA, sigma_V_sq_hat = NA, 
                num_active_weights = length(active_indices), 
                debug_info = "Frechet mean estimation failed"))
  }
  
  # 3. Calculate squared distances from active observations to the estimated Frechet mean (l_hat)
  d_sq_vals_active <- sapply(active_indices, function(i_act) {
    dist_fun_sq_calculator(Y_obj_list[[i_act]], l_hat, frechet_optns)
  })
  
  # 4. Estimate one-sided Frechet variance (V_hat)
  # V_hat = sum(w_i * d(Y_i, l_hat)^2) / sum(w_i)
  # Uses signed weights for the mean of squared distances (as per original implementation).
  sum_w_signed <- sum(s_in_weights_active)     
  V_hat        <- sum(s_in_weights_active * d_sq_vals_active, na.rm = TRUE) / sum_w_signed  
  
  # 5. Estimate variance of the Frechet variance estimator (sigma_V_sq_hat)
  # sigma_V_sq_hat = sum(|w_i| * (d_sq_i - V_hat)^2) / sum(|w_i|)
  # Uses absolute weights for the sum of squares part of variance calculation.
  sum_w_sq_abs <- sum(abs(s_in_weights_active)) 
  center       <- d_sq_vals_active - V_hat # Deviations of squared distances from V_hat
  
  sigma_V_sq_hat <- if (sum_w_sq_abs > 1e-12) { # Avoid division by zero if all weights are tiny
    sum(abs(s_in_weights_active) * center^2, na.rm = TRUE) / sum_w_sq_abs
  } else {
    NA # Indicate failure or insufficient data/weights
  }
  
  # Apply a floor to sigma_V_sq_hat to prevent issues with zero or negative variance in subsequent calculations
  if (is.na(sigma_V_sq_hat) || sigma_V_sq_hat < 1e-12) { 
    sigma_V_sq_hat <- 1e-12
  }
  
  return(list(l_hat = l_hat, V_hat = V_hat, sigma_V_sq_hat = sigma_V_sq_hat, 
              num_active_weights = length(active_indices), debug_info = NULL))
}

# Determines a heuristic range [min_bw, max_bw] for bandwidth selection
# for one side of an evaluation point 'c_val_target'.
# This is styled after bandwidth selection logic in the 'frechet' package.
.set_bw_range_frechet_style <- function(xin_side,          # Scalar covariate values on one side of c_val
                                        c_val_target,      # The evaluation point
                                        kernel_type,       # Character name of kernel (e.g., "gauss")
                                        X_scalar_global = NULL) { # All scalar covariate values (optional, for global SD fallback)
  
  # Handle cases with too few points on the side to reliably determine range
  if (length(xin_side) < 2) {
    warning("Too few points to determine CV bandwidth range robustly for one side. Using wide heuristic.")
    sd_global_val <- NA
    if (!is.null(X_scalar_global)) {
      sd_global_val <- stats::sd(X_scalar_global, na.rm = TRUE)
    }
    if (is.na(sd_global_val) || sd_global_val == 0) {
      sd_global_val <- 1 # Default SD if NA or zero to avoid issues
    }
    min_bw <- 0.01 * sd_global_val
    max_bw <- 0.5  * sd_global_val
    return(list(min = max(1e-5, min_bw),                # Ensure positive min_bw
                max = max(max_bw, min_bw * 1.5)))     # Ensure max is at least slightly larger than min
  }
  
  xinSt <- unique(sort(xin_side)) # Sorted unique points on the side
  
  # (The following xout_eff_min/max were in original but not used for bw range; kept for context if needed later)
  # xout_eff_min <- max(min(xinSt), c_val_target - 0.01 * diff(range(xinSt, na.rm = TRUE)))
  # xout_eff_max <- min(max(xinSt), c_val_target + 0.01 * diff(range(xinSt, na.rm = TRUE)))
  # if (xout_eff_min > xout_eff_max) { # Ensure min <= max
  #   xout_eff_min <- min(xinSt)
  #   xout_eff_max <- max(xinSt)
  # }
  
  # Minimum bandwidth heuristic based on minimum spacing between unique points
  min_spacing <- if (length(xinSt) > 1) min(diff(xinSt)) else 0.001 # Default if only one unique point (or <2 original)
  bw_min_heuristic <- min_spacing * 1.5
  if (is.na(bw_min_heuristic) || bw_min_heuristic == 0) {
    bw_min_heuristic <- 0.001 # Fallback
  }
  
  # Adjust for kernel type (some kernels have wider effective support for the same 'h')
  kernel_scale_factor <- (ifelse(tolower(kernel_type) == "gauss", 3, 1) *       # Gaussian kernel is wider
                            ifelse(tolower(kernel_type) == "gausvar", 2.5, 1))  # Variant Gaussian
  
  bw.min <- max(1e-5, bw_min_heuristic / kernel_scale_factor) # Ensure positive
  
  # Maximum bandwidth heuristic based on range or SD of points on the side
  range_xin_side <- diff(range(xin_side, na.rm = TRUE))
  sd_xin_side    <- stats::sd(xin_side, na.rm = TRUE)
  
  # Fallback for range_xin_side if it's zero or NA (e.g., all points are identical on this side)
  if (is.na(range_xin_side) || range_xin_side == 0) {
    if (!is.na(sd_xin_side) && sd_xin_side > 0) {
      range_xin_side <- sd_xin_side * 2 # Use SD if range is zero but SD is not
    } else {
      # Try global SD if side-specific SD is also problematic
      if (!is.null(X_scalar_global)) {
        sd_global_val_inner <- stats::sd(X_scalar_global, na.rm = TRUE)
        if (!is.na(sd_global_val_inner) && sd_global_val_inner > 0) {
          range_xin_side <- sd_global_val_inner * 2
        } else {
          range_xin_side <- 1 # Ultimate fallback if all else fails
        }
      } else {
        range_xin_side <- 1 # Ultimate fallback
      }
    }
  }
  if (is.na(range_xin_side) || range_xin_side == 0) {
    range_xin_side <- 1 # Ensure it's positive for bw.max calculation
  }
  
  bw.max <- max(bw.min * 1.5, range_xin_side / 3, na.rm = TRUE) # Ensure max is at least 1.5*min and related to data range
  
  # Ensure bw.max is strictly greater than bw.min and both are finite and reasonable
  if (bw.max < bw.min) {
    # This original condition `bw.min > bw.max * 1.5` seems unlikely if bw.max < bw.min.
    # The primary goal is to ensure bw.max > bw.min.
    if (bw.min > bw.max * 1.5 && is.finite(bw.min) && is.finite(bw.max)) { # Original condition
      bw.max <- bw.min * 1.01 # Make max slightly larger
    } else {
      bw.max <- (bw.min %||% 0.001) * 1.5 # Ensure max is at least 1.5 times min
    }
  }
  
  bw.min <- if (is.finite(bw.min)) bw.min else 1e-5 # Ensure min is finite
  bw.max <- if (is.finite(bw.max)) bw.max else bw.min * 2 # Ensure max is finite and larger than min
  
  if (bw.max <= bw.min) { # Final check to ensure max > min
    bw.max <- bw.min * 1.5 
  }
  
  return(list(min = bw.min, max = bw.max))
}

# Performs K-fold cross-validation to select a common bandwidth 'h' for estimating
# one-sided Frechet quantities on both sides of an evaluation point 'c_val'.
# The CV aims to minimize the sum of squared prediction errors on validation folds.
._cv_for_frechet_jump_h <- function(Y_obj_list,             # List of Frechet objects
                                    X_scalar,               # Scalar covariate
                                    c_val,                  # Evaluation point for the jump
                                    metric_space_type,      # Type of metric space
                                    kernel_frechet_char,    # Character name of kernel (e.g., "gauss")
                                    actual_kernel_function, # The kernel function itself
                                    frechet_optns,          # Options for Frechet mean estimation
                                    dist_fun_sq_calculator, # Function to calculate squared distance
                                    K_folds = 5,            # Number of CV folds
                                    n_bw_candidates = 10,   # Number of bandwidths to test
                                    min_eff_points_cv_fold = 3, # Min points in a training fold side
                                    min_bw_heuristic_factor = 0.001) { # Factor for min bandwidth heuristic
  
  n_total        <- length(X_scalar)
  X_scalar_range <- diff(range(X_scalar, na.rm = TRUE))
  if (!is.finite(X_scalar_range) || X_scalar_range == 0) {
    X_scalar_range <- 1 # Avoid division by zero or NA issues if range is degenerate
  }
  
  # Identify points to the left and right of c_val
  indices_left_all  <- which(X_scalar < c_val)
  indices_right_all <- which(X_scalar >= c_val)
  n_on_left_side    <- length(indices_left_all)
  n_on_right_side   <- length(indices_right_all)
  
  # Minimum number of points required on each side for CV to be feasible and somewhat robust
  # (K_folds-1)/K_folds of points are in training set. We need at least min_eff_points_cv_fold there.
  min_points_for_side_in_cv <- max(K_folds, # Need at least K_folds points to make K distinct folds
                                   ceiling(min_eff_points_cv_fold * K_folds / (K_folds - 1)), 
                                   5) # Absolute minimum
  
  # --- Fallback to Heuristic Bandwidth (if not enough data for CV) ---
  if (n_on_left_side < min_points_for_side_in_cv || n_on_right_side < min_points_for_side_in_cv) {
    warning(paste("Not enough points for robust CV on at least one side (left:", n_on_left_side,
                  ", right:", n_on_right_side, "; need >=", round(min_points_for_side_in_cv, 1),
                  " on each side). Using heuristic bandwidth for the jump."))
    
    X_around_c <- X_scalar[c(indices_left_all, indices_right_all)] # Data near c_val
    if (length(X_around_c) < 2) { # Not enough data even for simple bw.nrd0
      sd_X_global_fallback <- stats::sd(X_scalar, na.rm = TRUE)
      if (is.na(sd_X_global_fallback) || sd_X_global_fallback == 0) sd_X_global_fallback <- 1
      return(max(1e-5, 0.1 * sd_X_global_fallback, min_bw_heuristic_factor * X_scalar_range))
    }
    
    h_heuristic <- tryCatch(stats::bw.nrd0(X_around_c), error = function(e) NA) # Silverman's rule of thumb
    if (is.na(h_heuristic) || h_heuristic <= 0) { # If bw.nrd0 fails or gives non-positive
      h_heur_left <- NA; h_heur_right <- NA
      if (n_on_left_side >= 2)  h_heur_left  <- tryCatch(stats::bw.nrd0(X_scalar[indices_left_all]), error = function(e) NA)
      if (n_on_right_side >= 2) h_heur_right <- tryCatch(stats::bw.nrd0(X_scalar[indices_right_all]), error = function(e) NA)
      
      valid_heur <- c(h_heur_left, h_heur_right)
      valid_heur <- valid_heur[is.finite(valid_heur) & (valid_heur > 0)] # Filter out NA/non-positive
      
      if (length(valid_heur) > 0) {
        h_heuristic <- mean(valid_heur, na.rm = TRUE)
      } else { # If still no valid heuristic, use a fraction of global SD
        h_heuristic <- 0.1 * stats::sd(X_scalar, na.rm = TRUE)
        if (is.na(h_heuristic) || h_heuristic <= 0) h_heuristic <- 0.1 # Ultimate fallback value
      }
    }
    return(max(1e-5, h_heuristic, min_bw_heuristic_factor * X_scalar_range, na.rm = TRUE)) # Ensure positive
  }
  
  # --- Determine Bandwidth Candidate Range for CV ---
  bw_range_left <- .set_bw_range_frechet_style(xin_side        = X_scalar[indices_left_all],
                                               c_val_target    = c_val,
                                               kernel_type     = kernel_frechet_char,
                                               X_scalar_global = X_scalar)
  bw_range_right <- .set_bw_range_frechet_style(xin_side        = X_scalar[indices_right_all],
                                                c_val_target    = c_val,
                                                kernel_type     = kernel_frechet_char,
                                                X_scalar_global = X_scalar)
  
  # Combine ranges: take max of mins, min of maxs, bounded by heuristic factor of total range
  min_h_cand <- max(bw_range_left$min, bw_range_right$min, 
                    min_bw_heuristic_factor * X_scalar_range, na.rm = TRUE) 
  max_h_cand <- min(bw_range_left$max, bw_range_right$max, na.rm = TRUE)
  
  # Handle invalid or too narrow combined bandwidth range
  if (!is.finite(min_h_cand) || !is.finite(max_h_cand) || min_h_cand >= max_h_cand) {
    warning(paste("Combined CV: Invalid bandwidth range derived (min:",
                  signif(min_h_cand, 3), "max:", signif(max_h_cand, 3), 
                  "). Using global fallback range."))
    sd_X_global_fallback <- stats::sd(X_scalar, na.rm = TRUE)
    if (is.na(sd_X_global_fallback) || sd_X_global_fallback == 0) sd_X_global_fallback <- 1
    
    min_h_cand <- max(1e-5, 0.05 * sd_X_global_fallback, min_bw_heuristic_factor * X_scalar_range)
    max_h_cand <- max(min_h_cand * 1.5, 0.5 * sd_X_global_fallback) # Ensure max > min
  }
  
  bw_candidates <- seq(min_h_cand, max_h_cand, length.out = n_bw_candidates)
  bw_candidates <- unique(pmax(1e-5, bw_candidates)) # Ensure positive and unique candidates
  
  # --- K-Fold Cross-Validation Loop ---
  # Prepare K-fold assignments for left and right side data separately
  shuffled_indices_left_within_side  <- sample(1:n_on_left_side) # Shuffle indices *within* the left side subset
  folds_left                         <- cut(shuffled_indices_left_within_side, breaks = K_folds, labels = FALSE)
  shuffled_indices_right_within_side <- sample(1:n_on_right_side) # Shuffle indices *within* the right side subset
  folds_right                        <- cut(shuffled_indices_right_within_side, breaks = K_folds, labels = FALSE)
  
  # Calculate CV error for each candidate bandwidth
  cv_errors <- sapply(bw_candidates, function(h_cand) {
    fold_errors_h_sum          <- 0 # Sum of errors for this h_cand over all folds
    num_valid_folds_for_h_cand <- 0 # Count of folds that provided a valid error for this h_cand
    
    # Suppress warnings from inner estimations (e.g., too few points in a local window for locpoly) during CV
    suppressWarnings({ 
      for (k_fold_idx in 1:K_folds) {
        # Determine training and validation indices for the current fold (left side)
        validation_idx_left_in_shuffled <- which(folds_left == k_fold_idx)
        training_idx_left_in_shuffled   <- which(folds_left != k_fold_idx)
        # Map back to original indices in X_scalar and Y_obj_list
        original_indices_train_left_fold <- indices_left_all[shuffled_indices_left_within_side[training_idx_left_in_shuffled]]
        original_indices_val_left_fold   <- indices_left_all[shuffled_indices_left_within_side[validation_idx_left_in_shuffled]]
        
        # Determine training and validation indices for the current fold (right side)
        validation_idx_right_in_shuffled <- which(folds_right == k_fold_idx)
        training_idx_right_in_shuffled   <- which(folds_right != k_fold_idx)
        # Map back to original indices
        original_indices_train_right_fold <- indices_right_all[shuffled_indices_right_within_side[training_idx_right_in_shuffled]]
        original_indices_val_right_fold   <- indices_right_all[shuffled_indices_right_within_side[validation_idx_right_in_shuffled]]
        
        # Skip fold if not enough training data on either side for this h_cand
        if (length(original_indices_train_left_fold) < min_eff_points_cv_fold ||
            length(original_indices_train_right_fold) < min_eff_points_cv_fold) {
          next # Move to the next fold
        }
        
        # Estimate Frechet mean on the left side using training data of this fold
        l_hat_train_left_at_c <- NULL
        if (length(original_indices_train_left_fold) > 0) {
          s_in_train_weights_left <- .calculate_one_sided_locpoly_weights(
            x_obs = X_scalar[original_indices_train_left_fold], # Pass subset of X_scalar for training
            c_val = c_val, h_val = h_cand, K_fun = actual_kernel_function, side = "-", 
            n     = length(original_indices_train_left_fold) # n for local poly is size of this training subset
          )
          active_train_left_indices_in_subset <- which(s_in_train_weights_left != 0)
          
          if (length(active_train_left_indices_in_subset) >= min_eff_points_cv_fold) {
            l_hat_train_left_at_c <- .estimate_weighted_frechet_mean(
              Y_obj_list          = Y_obj_list[original_indices_train_left_fold], # Pass subset of Y_obj_list
              active_idx          = active_train_left_indices_in_subset,      # Indices relative to this Y_obj_list subset
              s_in_weights_active = s_in_train_weights_left[active_train_left_indices_in_subset],
              metric_space_type   = metric_space_type, 
              frechet_optns       = frechet_optns
            )
          }
          if (is.null(l_hat_train_left_at_c)) { next } # Skip if mean estimation failed for this fold/h_cand
        }
        
        # Estimate Frechet mean on the right side using training data of this fold
        l_hat_train_right_at_c <- NULL
        if (length(original_indices_train_right_fold) > 0) {
          s_in_train_weights_right <- .calculate_one_sided_locpoly_weights(
            x_obs = X_scalar[original_indices_train_right_fold], # Pass subset of X_scalar
            c_val = c_val, h_val = h_cand, K_fun = actual_kernel_function, side = "+", 
            n     = length(original_indices_train_right_fold)  # n for local poly is size of this training subset
          )
          active_train_right_indices_in_subset <- which(s_in_train_weights_right != 0)
          
          if (length(active_train_right_indices_in_subset) >= min_eff_points_cv_fold) {
            l_hat_train_right_at_c <- .estimate_weighted_frechet_mean(
              Y_obj_list          = Y_obj_list[original_indices_train_right_fold], # Pass subset of Y_obj_list
              active_idx          = active_train_right_indices_in_subset,       # Indices relative to this subset
              s_in_weights_active = s_in_train_weights_right[active_train_right_indices_in_subset],
              metric_space_type   = metric_space_type, 
              frechet_optns       = frechet_optns
            )
          }
          if (is.null(l_hat_train_right_at_c)) { next } # Skip if mean estimation failed
        }
        
        # If either mean estimation failed for the training fold, skip this fold for this h_cand
        if (is.null(l_hat_train_left_at_c) || is.null(l_hat_train_right_at_c)) {
          next 
        }
        
        # Calculate sum of squared prediction errors on validation sets for this fold
        current_fold_error_sum <- 0
        if (length(original_indices_val_left_fold) > 0) {
          errors_this_fold_val_left <- sapply(original_indices_val_left_fold, function(val_idx) {
            dist_fun_sq_calculator(Y_obj_list[[val_idx]], l_hat_train_left_at_c, frechet_optns)
          })
          current_fold_error_sum <- current_fold_error_sum + sum(errors_this_fold_val_left, na.rm = TRUE)
        }
        if (length(original_indices_val_right_fold) > 0) {
          errors_this_fold_val_right <- sapply(original_indices_val_right_fold, function(val_idx) {
            dist_fun_sq_calculator(Y_obj_list[[val_idx]], l_hat_train_right_at_c, frechet_optns)
          })
          current_fold_error_sum <- current_fold_error_sum + sum(errors_this_fold_val_right, na.rm = TRUE)
        }
        
        # Accumulate error if validation data was present and processed for this fold
        if (length(original_indices_val_left_fold) > 0 || length(original_indices_val_right_fold) > 0) {
          fold_errors_h_sum          <- fold_errors_h_sum + current_fold_error_sum
          num_valid_folds_for_h_cand <- num_valid_folds_for_h_cand + 1
        }
      } # End k_fold_idx loop (inner CV loop for folds)
    }) # End suppressWarnings
    
    if (num_valid_folds_for_h_cand == 0) return(Inf) # No valid folds produced an error for this h_cand
    return(fold_errors_h_sum / num_valid_folds_for_h_cand) # Average error for this h_cand across valid folds
  }) # End sapply over bw_candidates (outer CV loop for bandwidths)
  
  # --- Select Best Bandwidth and Handle CV Failure ---
  best_bw_idx <- which.min(cv_errors)
  
  # Handle cases where CV fails to find a good bandwidth (e.g., all errors are Inf)
  if (length(best_bw_idx) == 0 || !is.finite(cv_errors[best_bw_idx[1]]) || cv_errors[best_bw_idx[1]] == Inf) {
    warning_msg_cv_fail <- paste(
      "Combined K-Fold CV failed to find a valid minimum error. Min error found:",
      if (length(best_bw_idx) > 0 && length(cv_errors) >= best_bw_idx[1]) {
        signif(cv_errors[best_bw_idx[1]], 3) 
      } else { "NA" },
      ". Using heuristic fallback."
    )
    warning(warning_msg_cv_fail)
    
    # Re-apply heuristic fallback (similar to the one at the beginning of the function)
    X_around_c <- X_scalar[c(indices_left_all, indices_right_all)]
    if (length(X_around_c) < 2) {
      sd_X_global_fallback <- stats::sd(X_scalar, na.rm = TRUE)
      if (is.na(sd_X_global_fallback) || sd_X_global_fallback == 0) sd_X_global_fallback <- 1
      return(max(1e-5, 0.1 * sd_X_global_fallback, min_bw_heuristic_factor * X_scalar_range))
    }
    
    h_heuristic <- tryCatch(stats::bw.nrd0(X_around_c), error = function(e) NA)
    if (is.na(h_heuristic) || h_heuristic <= 0) {
      h_heur_left <- NA; h_heur_right <- NA
      if (n_on_left_side >= 2)  h_heur_left  <- tryCatch(stats::bw.nrd0(X_scalar[indices_left_all]), error = function(e) NA)
      if (n_on_right_side >= 2) h_heur_right <- tryCatch(stats::bw.nrd0(X_scalar[indices_right_all]), error = function(e) NA)
      
      valid_heur <- c(h_heur_left, h_heur_right)
      valid_heur <- valid_heur[is.finite(valid_heur) & (valid_heur > 0)]
      
      if (length(valid_heur) > 0) {
        h_heuristic <- mean(valid_heur, na.rm = TRUE)
      } else {
        h_heuristic <- 0.1 * stats::sd(X_scalar, na.rm = TRUE)
        if (is.na(h_heuristic) || h_heuristic <= 0) h_heuristic <- 0.1
      }
    }
    return(max(1e-5, h_heuristic, min_bw_heuristic_factor * X_scalar_range, na.rm = TRUE))
  }
  
  return(bw_candidates[best_bw_idx[1]]) # Return the best bandwidth found by CV
}