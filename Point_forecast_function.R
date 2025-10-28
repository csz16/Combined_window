# =====================================================================
# 0) SETUP
# =====================================================================
suppressPackageStartupMessages({
  library(ftsa)
  library(dplyr)
})

# Shorthands
ftsm          <- ftsa::ftsm
forecast.ftsm <- ftsa::forecast.ftsm

# =====================================================================
# 1) FORECAST GENERATORS (robust windowing)
# =====================================================================
rolling_forecast_iteration <- function(data_smooth, start_idx, end_idx, h, num_components) {
  train_window <- data_smooth[, start_idx:end_idx, drop = FALSE]
  model_fit <- tryCatch(ftsm(fts(1:nrow(data_smooth), train_window), order = num_components),
                        error = function(e) NULL)
  if (is.null(model_fit)) return(rep(NA_real_, nrow(data_smooth)))
  forecast_result <- tryCatch(forecast.ftsm(model_fit, h = h), error = function(e) NULL)
  if (is.null(forecast_result)) return(rep(NA_real_, nrow(data_smooth)))
  forecast_result$mean$y[, h, drop = TRUE]
}

expanding_forecast_iteration <- function(data_smooth, end_idx, h, num_components) {
  train_window <- data_smooth[, 1:end_idx, drop = FALSE]
  model_fit <- tryCatch(ftsm(fts(1:nrow(data_smooth), train_window), order = num_components),
                        error = function(e) NULL)
  if (is.null(model_fit)) return(rep(NA_real_, nrow(data_smooth)))
  forecast_result <- tryCatch(forecast.ftsm(model_fit, h = h), error = function(e) NULL)
  if (is.null(forecast_result)) return(rep(NA_real_, nrow(data_smooth)))
  forecast_result$mean$y[, h, drop = TRUE]
}

# =====================================================================
# 2) BETA SELECTION (POINT FORECAST; RMSE or MAE)
# =====================================================================
# TRUE  -> evaluate on RAW (exp forecasts, compare to validation_data_raw)
# FALSE -> evaluate on SMOOTH/log (compare to validation_data_smooth)
select_optimal_beta_point <- function(rolling_forecast_matrix_val,   # SMOOTH/log forecasts [age x n_val]
                                      expanding_forecast_matrix_val, # SMOOTH/log forecasts [age x n_val]
                                      validation_data_raw,           # RAW actuals [age x n_val]
                                      validation_data_smooth,        # SMOOTH/log actuals [age x n_val]
                                      apply_exp_transform,
                                      objective = c("rmse","mae"),
                                      grid_points = 2001,
                                      beta_diagnostics = FALSE) {
  objective <- match.arg(objective)
  
  # Route evaluation scale
  if (isTRUE(apply_exp_transform)) {
    roll_eval <- exp(rolling_forecast_matrix_val)
    expd_eval <- exp(expanding_forecast_matrix_val)
    actuals   <- validation_data_raw
  } else {
    roll_eval <- rolling_forecast_matrix_val
    expd_eval <- expanding_forecast_matrix_val
    actuals   <- validation_data_smooth
  }
  
  valid <- is.finite(roll_eval) & is.finite(expd_eval) & is.finite(actuals)
  if (!any(valid)) {
    warning("No valid cells for validation objective; defaulting beta=0.5.")
    return(list(beta = 0.5, objective_value = NA_real_, objective_name = objective,
                which_method = NA_integer_, at_boundary = NA, diagnostics_curve = NULL))
  }
  
  # Objective functions
  mse_of_beta <- function(b) {
    pred <- b * roll_eval + (1 - b) * expd_eval
    v <- (actuals - pred)^2
    mean(v[valid])
  }
  mae_of_beta <- function(b) {
    pred <- b * roll_eval + (1 - b) * expd_eval
    v <- abs(actuals - pred)
    mean(v[valid])
  }
  
  if (objective == "rmse") {
    o <- optimise(f = mse_of_beta, interval = c(0, 1))
    best_beta <- o$minimum
    best_val  <- sqrt(o$objective)  # report RMSE
    which_m   <- 1L
  } else {
    # MAE is piecewise-linear; use dense grid search
    grid <- seq(0, 1, length.out = grid_points)
    vals <- vapply(grid, mae_of_beta, numeric(1))
    k <- which.min(vals)
    best_beta <- grid[k]
    best_val  <- vals[k]
    which_m   <- 2L
  }
  
  at_boundary <- (best_beta <= 1e-6) || (best_beta >= 1 - 1e-6)
  
  diag_df <- NULL
  if (isTRUE(beta_diagnostics)) {
    grid <- seq(0, 1, length.out = grid_points)
    obj_vals <- if (objective == "rmse") sqrt(vapply(grid, mse_of_beta, numeric(1)))
    else vapply(grid, mae_of_beta, numeric(1))
    diag_df <- data.frame(beta = grid, objective = obj_vals)
  }
  
  list(
    beta = as.numeric(best_beta),
    objective_value = best_val,
    objective_name  = objective,
    which_method    = which_m,
    at_boundary     = at_boundary,
    diagnostics_curve = diag_df
  )
}

# =====================================================================
# 3) TEST EVALUATION (point-only)
# =====================================================================
compare_test_performance_point <- function(forecast_list, actual_test_data, apply_exp_transform) {
  method_names <- names(forecast_list)
  results <- data.frame()
  scale_msg <- if (apply_exp_transform) "RAW (exp/additive)" else "SMOOTH/log (additive)"
  cat(sprintf("Comparing performance on test set (%s)...\n", scale_msg))
  
  for (method in method_names) {
    pf <- forecast_list[[method]]
    valid_idx <- is.finite(actual_test_data) & is.finite(pf$mean)
    actuals <- as.vector(actual_test_data[valid_idx])
    preds   <- as.vector(pf$mean[valid_idx])
    
    rmse <- if (length(actuals) > 0) sqrt(mean((actuals - preds)^2)) else NA_real_
    mae  <- if (length(actuals) > 0) mean(abs(actuals - preds)) else NA_real_
    
    results <- rbind(results, data.frame(
      Method = method, RMSE = rmse, MAE = mae,
      Exp_Transform = apply_exp_transform,
      row.names = NULL
    ))
  }
  results
}

# =====================================================================
# 4) MAIN WRAPPER (POINT FORECAST ONLY)
# =====================================================================
forecast_with_point_beta_optimization <- function(data_smooth,
                                                  data_raw,
                                                  end_train_index,
                                                  h,
                                                  num_components,
                                                  apply_exp_transform = TRUE,  # TRUE->RAW; FALSE->SMOOTH/log
                                                  seed = 42,
                                                  objective = c("rmse","mae"),
                                                  beta_diagnostics = FALSE,
                                                  grid_points = 2001) {
  objective <- match.arg(objective)
  
  set.seed(seed)
  
  num_rows    <- nrow(data_smooth)
  num_periods <- ncol(data_smooth)
  train_len   <- end_train_index
  
  # All potentially-available targets for this h
  available_all <- (end_train_index + h):num_periods
  if (length(available_all) < 2L) stop("Not enough data for validation and test sets.")
  
  # Baseline total targets at h = 1  (e.g., 40 when end_train_index = ncol - 40)
  total_start <- num_periods - end_train_index
  base_pairs  <- floor(total_start / 2L)          # e.g., 40 -> 20
  
  # Pairs per side follow 20,19,18,...,1 for h = 1.. when total_start permits
  pairs_h <- max(base_pairs - (h - 1L), 1L)
  val_len  <- pairs_h
  test_len <- pairs_h
  
  # Positions inside available_all:
  val_pos <- seq_len(val_len)
  test_start_pos <- val_len + h
  test_pos <- seq.int(test_start_pos, test_start_pos + test_len - 1L)
  
  if (max(test_pos) > length(available_all)) {
    stop(sprintf("Split overflow: need %d positions, but only %d available (h=%d).",
                 max(test_pos), length(available_all), h))
  }
  
  validation_target_indices <- available_all[val_pos]
  test_target_indices       <- available_all[test_pos]
  
  cat(sprintf(
    "Split (h=%d): total_start=%d | available=%d | val=%d (%d..%d) | test=%d (%d..%d) | gap=%d\n",
    h, total_start, length(available_all),
    val_len,  min(validation_target_indices), max(validation_target_indices),
    test_len, min(test_target_indices),       max(test_target_indices),
    h - 1L
  ))
  
  # ===== Validation matrices =====
  validation_data_raw_actuals    <- data_raw[,    validation_target_indices,    drop = FALSE]
  validation_data_smooth_actuals <- data_smooth[, validation_target_indices,    drop = FALSE]
  
  # --- Validation forecasts (SMOOTH/log)
  cat("\n--- Starting Validation Phase (point forecast) ---\n")
  forecast_rates_rolling_val   <- matrix(NA_real_, nrow = num_rows, ncol = val_len)
  forecast_rates_expanding_val <- matrix(NA_real_, nrow = num_rows, ncol = val_len)
  for (i in 1:val_len) {
    end_idx <- validation_target_indices[i] - h
    start_idx <- max(1, end_idx - train_len + 1)
    forecast_rates_rolling_val[, i]   <- rolling_forecast_iteration(data_smooth, start_idx, end_idx, h, num_components)
    forecast_rates_expanding_val[, i] <- expanding_forecast_iteration(data_smooth, end_idx, h, num_components)
  }
  
  # --- Beta selection (point forecast objective)
  cat("\n--- Optimizing beta (", objective, ") ---\n", sep = "")
  res_optimal <- select_optimal_beta_point(
    rolling_forecast_matrix_val   = forecast_rates_rolling_val,
    expanding_forecast_matrix_val = forecast_rates_expanding_val,
    validation_data_raw           = validation_data_raw_actuals,
    validation_data_smooth        = validation_data_smooth_actuals,
    apply_exp_transform           = apply_exp_transform,
    objective                     = objective,
    grid_points                   = grid_points,
    beta_diagnostics              = beta_diagnostics
  )
  optimal_beta <- res_optimal$beta
  cat(sprintf("Optimal beta = %.4f%s | %s = %.6f\n",
              optimal_beta,
              if (isTRUE(res_optimal$at_boundary)) " (boundary)" else "",
              toupper(objective), res_optimal$objective_value))
  
  # --- Test forecasts (SMOOTH/log)
  cat("\n--- Building test point forecasts on evaluation scale ---\n")
  test_data_raw_actuals    <- data_raw[,    test_target_indices, drop = FALSE]
  test_data_smooth_actuals <- data_smooth[, test_target_indices, drop = FALSE]
  
  n_test <- length(test_target_indices)
  fc_rolling_test   <- matrix(NA_real_, nrow = num_rows, ncol = n_test)
  fc_expanding_test <- matrix(NA_real_, nrow = num_rows, ncol = n_test)
  for (i in 1:n_test) {
    end_idx <- test_target_indices[i] - h
    start_idx <- max(1, end_idx - train_len + 1)
    fc_rolling_test[, i]   <- rolling_forecast_iteration(data_smooth, start_idx, end_idx, h, num_components)
    fc_expanding_test[, i] <- expanding_forecast_iteration(data_smooth, end_idx, h, num_components)
  }
  
  # ---------- Build point means on the evaluation scale ----------
  if (isTRUE(apply_exp_transform)) {
    roll_mean_raw <- exp(fc_rolling_test)
    expd_mean_raw <- exp(fc_expanding_test)
    
    combined_opt_mean_raw <- optimal_beta * roll_mean_raw + (1 - optimal_beta) * expd_mean_raw
    fix_beta <- 0.5
    combined_fix_mean_raw <- fix_beta * roll_mean_raw + (1 - fix_beta) * expd_mean_raw
    
    test_forecast_list_for_eval <- list(
      rolling          = list(mean = roll_mean_raw),
      expanding        = list(mean = expd_mean_raw),
      combined_optimal = list(mean = combined_opt_mean_raw),
      fix_beta         = list(mean = combined_fix_mean_raw)
    )
    actuals_for_eval <- test_data_raw_actuals
    
    # Also provide RAW means (same object, for convenience)
    test_forecast_list_raw <- test_forecast_list_for_eval
    
  } else {
    # SMOOTH/log scale evaluation
    roll_mean_sm <- fc_rolling_test
    expd_mean_sm <- fc_expanding_test
    
    combined_opt_mean_sm <- optimal_beta * roll_mean_sm + (1 - optimal_beta) * expd_mean_sm
    fix_beta <- 0.5
    combined_fix_mean_sm <- fix_beta * roll_mean_sm + (1 - fix_beta) * expd_mean_sm
    
    test_forecast_list_for_eval <- list(
      rolling          = list(mean = roll_mean_sm),
      expanding        = list(mean = expd_mean_sm),
      combined_optimal = list(mean = combined_opt_mean_sm),
      fix_beta         = list(mean = combined_fix_mean_sm)
    )
    actuals_for_eval <- test_data_smooth_actuals
    
    # For convenience, also provide RAW (exp) versions
    test_forecast_list_raw <- lapply(test_forecast_list_for_eval, function(m) {
      list(mean = exp(m$mean))
    })
  }
  
  cat(sprintf("\nEvaluating on %s scale...\n",
              if (apply_exp_transform) "RAW (exp/additive)" else "SMOOTH/log (additive)"))
  
  test_results <- compare_test_performance_point(
    forecast_list       = test_forecast_list_for_eval,
    actual_test_data    = actuals_for_eval,
    apply_exp_transform = apply_exp_transform
  )
  
  cat("\n--- Function Finished (POINT FORECAST) ---\n")
  list(
    settings = list(h = h,
                    exp_transform_applied = apply_exp_transform),
    optimization_results = list(optimal_beta_model = res_optimal),
    test_phase_performance = test_results,
    # Always include RAW test means for convenience
    test_phase_forecasts = c(list(actuals_raw = test_data_raw_actuals,
                                  actuals_smooth = test_data_smooth_actuals),
                             test_forecast_list_raw)
  )
}

# =====================================================================
# 5) EXAMPLE USAGE (POINT FORECAST)
# =====================================================================
# Example objects via Female mortality from AUS when h = 1
smooth_mort_f <- demog_list$AUS$Female_smooth$rate[[1]]
mort_f        <- demog_list$AUS$Female_raw$rate[[1]]

final_results <- forecast_with_point_beta_optimization(
  data_smooth         = log(smooth_mort_f),  # log SMOOTH matrix
  data_raw            = log(mort_f),              # log RAW matrix
  end_train_index     = 61,
  h                   = 1,
  num_components      = 6,
  apply_exp_transform = FALSE,                # FALSE => evaluate on SMOOTH/log
  seed                = 42,
  objective           = "mae"               # or "mae"
)

