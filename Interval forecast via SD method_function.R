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
# 2) METRICS
# =====================================================================
# Interval Score (TEST ONLY) — metrics fixed: ECP, CPD
interval_score <- function(holdout, lb, ub, alpha) {
  holdout <- as.vector(holdout); lb <- as.vector(lb); ub <- as.vector(ub)
  valid <- is.finite(holdout) & is.finite(lb) & is.finite(ub)
  if (!any(valid)) return(c(ECP = NA, CPD = NA, mean_score = Inf))
  
  holdout <- holdout[valid]; lb <- lb[valid]; ub <- ub[valid]
  below <- holdout < lb; above <- holdout > ub
  
  # ECP: empirical coverage; nominal coverage is (1 - alpha)
  ecp <- 1 - (sum(below) + sum(above)) / length(holdout)
  cpd <- abs(ecp - (1 - alpha))
  
  score <- (ub - lb) + (2/alpha) * ((lb - holdout) * below + (holdout - ub) * above)
  c(ECP = ecp, CPD = cpd, mean_score = mean(score))
}

# SD across age, with fallback
calculate_fts_sd <- function(residual_matrix, sd_fts_type = "coordinate") {
  if (any(!is.finite(residual_matrix))) residual_matrix[!is.finite(residual_matrix)] <- NA
  resi_fts <- fts(x = 0:(nrow(residual_matrix) - 1), y = residual_matrix)
  sd_vector <- tryCatch(
    ftsa::sd.fts(resi_fts, method = sd_fts_type)$y,
    error = function(e) {
      warning("ftsa::sd.fts failed: ", conditionMessage(e), ". Falling back to row SDs (na.rm=TRUE).")
      apply(residual_matrix, 1, sd, na.rm = TRUE)
    }
  )
  sd_vector[is.na(sd_vector)] <- 0
  sd_vector
}

# =====================================================================
# 3) TUNING (POINTWISE CPD target; returns abs(ECP - alpha_level))
# =====================================================================
tune_para_objective_fn <- function(par, forecast_matrix, validation_data, sd_vector, alpha_level) {
  tune_para <- as.numeric(par[1])
  if (!is.finite(tune_para) || tune_para < 0) return(Inf)
  
  residual_matrix <- validation_data - forecast_matrix
  n_age <- nrow(residual_matrix); n_col <- ncol(residual_matrix)
  
  # ensure sd_vector length & finiteness
  sd_vector <- as.numeric(sd_vector)
  if (length(sd_vector) != n_age) sd_vector <- rep_len(sd_vector, n_age)
  sd_vector[!is.finite(sd_vector)] <- 0
  
  # vectorised inclusion check (additive bands)
  lower_mat <- matrix(-tune_para * sd_vector, nrow = n_age, ncol = n_col)
  upper_mat <- matrix( tune_para * sd_vector, nrow = n_age, ncol = n_col)
  
  valid <- is.finite(residual_matrix) & is.finite(lower_mat) & is.finite(upper_mat)
  if (!any(valid)) return(Inf)
  
  ind <- (residual_matrix >= lower_mat) & (residual_matrix <= upper_mat)
  ecp <- sum(ind[valid], na.rm = TRUE) / sum(valid)
  
  abs(ecp - alpha_level)  # CPD
}

# Route scale EXACTLY per rule:
# TRUE  -> EXP forecasts, compare to RAW (validation_data_raw)
# FALSE -> keep SMOOTH/log forecasts, compare to SMOOTH/log (validation_data_smooth)
optimize_tuning_parameter <- function(forecast_matrix_val,   # SMOOTH/log forecasts
                                      validation_data_raw,
                                      validation_data_smooth,
                                      alpha_level,                 # target coverage (e.g., 0.80)
                                      apply_exp_transform = TRUE,  # TRUE->RAW path, FALSE->SMOOTH path
                                      sd_fts_type = "coordinate",
                                      seed = 42) {
  # route the scale for tuning
  if (isTRUE(apply_exp_transform)) {
    forecasts_for_opt <- exp(forecast_matrix_val)   # RAW
    actuals_for_opt   <- validation_data_raw        # RAW
  } else {
    forecasts_for_opt <- forecast_matrix_val        # SMOOTH/log
    actuals_for_opt   <- validation_data_smooth     # SMOOTH/log
  }
  
  # residuals and SD across "age" via ftsa::sd.fts (EXACT), with safe fallback
  residuals <- actuals_for_opt - forecasts_for_opt
  sd_vector <- tryCatch({
    resi_fts <- fts(seq_len(nrow(residuals)), residuals)  # row index as age
    ftsa::sd.fts(resi_fts, method = sd_fts_type)$y
  }, error = function(e) {
    apply(residuals, 1, sd, na.rm = TRUE)
  })
  sd_vector[!is.finite(sd_vector)] <- 0
  
  set.seed(seed)
  dum_1 <- suppressWarnings(optim(par = 2, fn = tune_para_objective_fn, method = "Nelder-Mead",
                                  forecast_matrix = forecasts_for_opt,
                                  validation_data = actuals_for_opt,
                                  sd_vector = sd_vector, alpha_level = alpha_level))
  
  dum_2 <- optim(par = 2, fn = tune_para_objective_fn, method = "Brent",
                 lower = 0, upper = 1e3,
                 forecast_matrix = forecasts_for_opt,
                 validation_data = actuals_for_opt,
                 sd_vector = sd_vector, alpha_level = alpha_level)
  
  dum_3 <- optim(par = 2, fn = tune_para_objective_fn, method = "L-BFGS-B",
                 lower = 0, upper = 1e3,
                 forecast_matrix = forecasts_for_opt,
                 validation_data = actuals_for_opt,
                 sd_vector = sd_vector, alpha_level = alpha_level)
  
  dum_4 <- optimise(f = tune_para_objective_fn, interval = c(0, 1e3),
                    forecast_matrix = forecasts_for_opt,
                    validation_data = actuals_for_opt,
                    sd_vector = sd_vector, alpha_level = alpha_level)
  
  vals <- c(dum_1$value, dum_2$value, dum_3$value, dum_4$objective)
  k <- which.min(vals)
  best_par <- if (k == 1) as.numeric(dum_1$par) else
    if (k == 2) as.numeric(dum_2$par) else
      if (k == 3) as.numeric(dum_3$par) else
        as.numeric(dum_4$minimum)
  
  list(tune_para = best_par, sd_vector = sd_vector)
}

# =====================================================================
# 4) BETA SELECTION (four-method sweep; CPD or Interval Score) — FIXED
# =====================================================================
# TRUE  -> EXP forecasts, compare to RAW
# FALSE -> SMOOTH/log forecasts, compare to SMOOTH/log
select_optimal_params <- function(rolling_forecast_matrix_val,   # SMOOTH/log forecasts [age x n_val]
                                  expanding_forecast_matrix_val, # SMOOTH/log forecasts [age x n_val]
                                  validation_data_raw,           # RAW actuals [age x n_val]
                                  validation_data_smooth,        # SMOOTH/log actuals [age x n_val]
                                  alpha_level,                   # coverage target (e.g., 0.80)
                                  apply_exp_transform,
                                  tune_para_rolling, sd_vector_rolling,
                                  tune_para_expanding, sd_vector_expanding,
                                  objective = c("cpd","interval_score"),
                                  seed = 42,
                                  beta_diagnostics = FALSE,
                                  grid_points = 201,
                                  tie_breaker_eps = 0) {
  objective <- match.arg(objective)
  
  # Route scale
  if (isTRUE(apply_exp_transform)) {
    rolling_fc_opt   <- exp(rolling_forecast_matrix_val)
    expanding_fc_opt <- exp(expanding_forecast_matrix_val)
    actuals_for_opt  <- validation_data_raw
  } else {
    rolling_fc_opt   <- rolling_forecast_matrix_val
    expanding_fc_opt <- expanding_forecast_matrix_val
    actuals_for_opt  <- validation_data_smooth
  }
  
  n_age <- nrow(actuals_for_opt); n_col <- ncol(actuals_for_opt)
  
  # Clean sd vectors to length n_age
  sd_vector_rolling   <- rep_len(as.numeric(sd_vector_rolling),   n_age)
  sd_vector_expanding <- rep_len(as.numeric(sd_vector_expanding), n_age)
  sd_vector_rolling[!is.finite(sd_vector_rolling)]     <- 0
  sd_vector_expanding[!is.finite(sd_vector_expanding)] <- 0
  
  pi_alpha <- 1 - alpha_level  # miscoverage for interval-score term
  
  make_bounds <- function(beta) {
    mean_fc   <- beta * rolling_fc_opt + (1 - beta) * expanding_fc_opt
    width_vec <- beta * tune_para_rolling * sd_vector_rolling +
      (1 - beta) * tune_para_expanding * sd_vector_expanding
    list(
      mean  = mean_fc,
      lower = sweep(mean_fc, 1, width_vec, "-"),
      upper = sweep(mean_fc, 1, width_vec, "+"),
      width_vec = width_vec
    )
  }
  
  # CPD objective: |ECP - alpha_level|
  cpd_objective <- function(beta) {
    b <- make_bounds(beta)
    x <- actuals_for_opt; l <- b$lower; u <- b$upper
    valid <- is.finite(x) & is.finite(l) & is.finite(u)
    if (!any(valid)) return(Inf)
    covered <- (x >= l) & (x <= u)
    ecp <- sum(covered[valid]) / sum(valid)  # denominator uses valid cells only
    base <- abs(ecp - alpha_level)           # CPD
    if (tie_breaker_eps > 0) {
      w <- mean(b$width_vec[is.finite(b$width_vec)])
      base + tie_breaker_eps * w  # tiny width penalty to break flats
    } else base
  }
  
  # Interval score objective
  ints_objective <- function(beta) {
    b <- make_bounds(beta)
    x <- actuals_for_opt; l <- b$lower; u <- b$upper
    valid <- is.finite(x) & is.finite(l) & is.finite(u)
    if (!any(valid)) return(Inf)
    base_width <- u - l
    below <- (x < l); above <- (x > u)
    penalty <- (l - x) * below + (x - u) * above
    mean((base_width + (2 / pi_alpha) * penalty)[valid])
  }
  
  obj_fun <- if (objective == "cpd") cpd_objective else ints_objective
  
  set.seed(seed)
  dum_1 <- optimise(f = obj_fun, interval = c(0, 1))
  dum_2 <- optim(par = 0.5, fn = obj_fun, method = "Brent",    lower = 0, upper = 1)
  dum_3 <- optim(par = 0.5, fn = obj_fun, method = "L-BFGS-B", lower = 0, upper = 1)
  
  vals <- c(dum_1$objective, dum_2$value, dum_3$value)
  k <- which.min(vals)
  best_beta <- c(dum_1$minimum, as.numeric(dum_2$par), as.numeric(dum_3$par))[k]
  best_val  <- min(vals)
  at_boundary <- (best_beta <= 1e-6) || (best_beta >= 1 - 1e-6)
  
  diag_df <- NULL
  if (isTRUE(beta_diagnostics)) {
    grid <- seq(0, 1, length.out = grid_points)
    obj_vals <- vapply(grid, obj_fun, numeric(1))
    cov_vals <- vapply(grid, function(b) {
      bb <- make_bounds(b)
      x <- actuals_for_opt; l <- bb$lower; u <- bb$upper
      valid <- is.finite(x) & is.finite(l) & is.finite(u)
      if (!any(valid)) return(NA_real_)
      covered <- (x >= l) & (x <= u)
      sum(covered[valid]) / sum(valid)
    }, numeric(1))
    diag_df <- data.frame(beta = grid, objective = obj_vals, ecp = cov_vals)
  }
  
  list(
    beta = as.numeric(best_beta),
    tune_para_rolling   = tune_para_rolling,
    tune_para_expanding = tune_para_expanding,
    sd_vector_rolling   = sd_vector_rolling,
    sd_vector_expanding = sd_vector_expanding,
    objective_value     = best_val,
    objective_name      = objective,  # "cpd" or "interval_score"
    which_method        = k,
    at_boundary         = at_boundary,
    diagnostics_curve   = diag_df
  )
}

# =====================================================================
# 5) TEST EVALUATION (on evaluation scale)
# =====================================================================
compare_test_performance <- function(forecast_list, actual_test_data, pi_alpha, apply_exp_transform) {
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
    
    ecp <- NA_real_; cpd <- NA_real_; mean_score <- NA_real_
    if (!is.null(pf$lower) && !is.null(pf$upper)) {
      m <- interval_score(holdout = actual_test_data, lb = pf$lower, ub = pf$upper, alpha = pi_alpha)
      ecp        <- unname(m["ECP"])
      cpd        <- unname(m["CPD"])
      mean_score <- unname(m["mean_score"])
    }
    
    results <- rbind(results, data.frame(
      Method = method, RMSE = rmse, MAE = mae,
      ECP = ecp, CPD = cpd,
      Mean_Interval_Score = mean_score,
      Target_Coverage = 1 - pi_alpha,
      Alpha = pi_alpha,
      Exp_Transform = apply_exp_transform,
      row.names = NULL
    ))
  }
  results
}

# =====================================================================
# 6) MAIN WRAPPER (routing EXACTLY per rule)
# =====================================================================
forecast_with_interval_optimization <- function(data_smooth,
                                                data_raw,
                                                end_train_index,
                                                h,
                                                num_components,
                                                pi_level = 0.80,             # coverage target
                                                apply_exp_transform = TRUE,  # TRUE->RAW; FALSE->SMOOTH/log
                                                sd_fts_type = "coordinate",
                                                seed = 42,
                                                objective = c("interval_score","cpd"),
                                                beta_diagnostics = FALSE,
                                                tie_breaker_eps = 0,
                                                grid_points = 201) {
  objective <- match.arg(objective)
  pi_alpha <- 1 - pi_level
  
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
  cat("\n--- Starting Validation Phase ---\n")
  forecast_rates_rolling_val   <- matrix(NA_real_, nrow = num_rows, ncol = val_len)
  forecast_rates_expanding_val <- matrix(NA_real_, nrow = num_rows, ncol = val_len)
  for (i in 1:val_len) {
    end_idx <- validation_target_indices[i] - h
    start_idx <- max(1, end_idx - train_len + 1)
    forecast_rates_rolling_val[, i]   <- rolling_forecast_iteration(data_smooth, start_idx, end_idx, h, num_components)
    forecast_rates_expanding_val[, i] <- expanding_forecast_iteration(data_smooth, end_idx, h, num_components)
  }
  
  # --- Tune rolling & expanding (POINTWISE CPD target; routing per flag)
  cat("\n--- Optimizing tuning parameters (pointwise CPD target) ---\n")
  res_roll <- optimize_tuning_parameter(
    forecast_matrix_val     = forecast_rates_rolling_val,
    validation_data_raw     = validation_data_raw_actuals,
    validation_data_smooth  = validation_data_smooth_actuals,
    alpha_level             = pi_level,
    apply_exp_transform     = apply_exp_transform,
    sd_fts_type             = sd_fts_type,
    seed                    = seed
  )
  res_expand <- optimize_tuning_parameter(
    forecast_matrix_val     = forecast_rates_expanding_val,
    validation_data_raw     = validation_data_raw_actuals,
    validation_data_smooth  = validation_data_smooth_actuals,
    alpha_level             = pi_level,
    apply_exp_transform     = apply_exp_transform,
    sd_fts_type             = sd_fts_type,
    seed                    = seed
  )
  cat(sprintf("Rolling tune_para = %.4f | Expanding tune_para = %.4f\n",
              res_roll$tune_para, res_expand$tune_para))
  
  # --- Beta selection (uses same scale as tuning)
  cat("\n--- Optimizing beta (", objective, ") ---\n", sep = "")
  res_optimal <- select_optimal_params(
    rolling_forecast_matrix_val   = forecast_rates_rolling_val,
    expanding_forecast_matrix_val = forecast_rates_expanding_val,
    validation_data_raw           = validation_data_raw_actuals,
    validation_data_smooth        = validation_data_smooth_actuals,
    alpha_level                   = pi_level,
    apply_exp_transform           = apply_exp_transform,
    tune_para_rolling             = res_roll$tune_para,
    sd_vector_rolling             = res_roll$sd_vector,
    tune_para_expanding           = res_expand$tune_para,
    sd_vector_expanding           = res_expand$sd_vector,
    objective                     = objective,
    seed                          = seed,
    beta_diagnostics              = beta_diagnostics,
    grid_points                   = grid_points,
    tie_breaker_eps               = tie_breaker_eps
  )
  optimal_beta <- res_optimal$beta
  cat(sprintf("Optimal beta = %.4f%s\n",
              optimal_beta,
              if (isTRUE(res_optimal$at_boundary)) " (boundary)" else ""))
  
  # --- Test forecasts (SMOOTH/log)
  cat("\n--- Building test forecasts & intervals on evaluation scale ---\n")
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
  
  # ---------- Build means & intervals on the evaluation scale ----------
  if (isTRUE(apply_exp_transform)) {
    # Means on RAW scale
    roll_mean_raw <- exp(fc_rolling_test)
    expd_mean_raw <- exp(fc_expanding_test)
    
    # Widths tuned on RAW scale
    width_roll_raw <- res_roll$tune_para   * res_roll$sd_vector     # [age]
    width_expd_raw <- res_expand$tune_para * res_expand$sd_vector   # [age]
    
    # Combine means and widths
    combined_opt_mean_raw <- optimal_beta * roll_mean_raw + (1 - optimal_beta) * expd_mean_raw
    width_opt_raw <- optimal_beta * width_roll_raw + (1 - optimal_beta) * width_expd_raw
    
    fix_beta <- 0.5
    combined_fix_mean_raw <- fix_beta * roll_mean_raw + (1 - fix_beta) * expd_mean_raw
    width_fix_raw <- fix_beta * width_roll_raw + (1 - fix_beta) * width_expd_raw
    
    # Forecast list for evaluation (RAW, additive)
    test_forecast_list_for_eval <- list(
      rolling = list(
        mean  = roll_mean_raw,
        lower = sweep(roll_mean_raw, 1, width_roll_raw, "-"),
        upper = sweep(roll_mean_raw, 1, width_roll_raw, "+")
      ),
      expanding = list(
        mean  = expd_mean_raw,
        lower = sweep(expd_mean_raw, 1, width_expd_raw, "-"),
        upper = sweep(expd_mean_raw, 1, width_expd_raw, "+")
      ),
      combined_optimal = list(
        mean  = combined_opt_mean_raw,
        lower = sweep(combined_opt_mean_raw, 1, width_opt_raw, "-"),
        upper = sweep(combined_opt_mean_raw, 1, width_opt_raw, "+")
      ),
      fix_beta = list(
        mean  = combined_fix_mean_raw,
        lower = sweep(combined_fix_mean_raw, 1, width_fix_raw, "-"),
        upper = sweep(combined_fix_mean_raw, 1, width_fix_raw, "+")
      )
    )
    actuals_for_eval <- test_data_raw_actuals
    
    # Also define RAW forecasts for return convenience
    test_forecast_list_raw <- test_forecast_list_for_eval
    
  } else {
    # Means on SMOOTH/log scale
    roll_mean_sm <- fc_rolling_test
    expd_mean_sm <- fc_expanding_test
    
    # Widths tuned on SMOOTH/log scale
    width_roll_sm <- res_roll$tune_para   * res_roll$sd_vector
    width_expd_sm <- res_expand$tune_para * res_expand$sd_vector
    
    combined_opt_mean_sm <- optimal_beta * roll_mean_sm + (1 - optimal_beta) * expd_mean_sm
    width_opt_sm <- optimal_beta * width_roll_sm + (1 - optimal_beta) * width_expd_sm
    
    fix_beta <- 0.5
    combined_fix_mean_sm <- fix_beta * roll_mean_sm + (1 - fix_beta) * expd_mean_sm
    width_fix_sm <- fix_beta * width_roll_sm + (1 - fix_beta) * width_expd_sm
    
    test_forecast_list_for_eval <- list(
      rolling = list(
        mean  = roll_mean_sm,
        lower = sweep(roll_mean_sm, 1, width_roll_sm, "-"),
        upper = sweep(roll_mean_sm, 1, width_roll_sm, "+")
      ),
      expanding = list(
        mean  = expd_mean_sm,
        lower = sweep(expd_mean_sm, 1, width_expd_sm, "-"),
        upper = sweep(expd_mean_sm, 1, width_expd_sm, "+")
      ),
      combined_optimal = list(
        mean  = combined_opt_mean_sm,
        lower = sweep(combined_opt_mean_sm, 1, width_opt_sm, "-"),
        upper = sweep(combined_opt_mean_sm, 1, width_opt_sm, "+")
      ),
      fix_beta = list(
        mean  = combined_fix_mean_sm,
        lower = sweep(combined_fix_mean_sm, 1, width_fix_sm, "-"),
        upper = sweep(combined_fix_mean_sm, 1, width_fix_sm, "+")
      )
    )
    actuals_for_eval <- test_data_smooth_actuals
    
    # For convenience, also produce RAW (exp) versions for return
    test_forecast_list_raw <- lapply(test_forecast_list_for_eval, function(m) {
      list(mean = exp(m$mean), lower = exp(m$lower), upper = exp(m$upper))
    })
  }
  
  cat(sprintf("\nEvaluating on %s scale...\n",
              if (apply_exp_transform) "RAW (exp/additive)" else "SMOOTH/log (additive)"))
  
  test_results <- compare_test_performance(
    forecast_list       = test_forecast_list_for_eval,
    actual_test_data    = actuals_for_eval,
    pi_alpha            = pi_alpha,
    apply_exp_transform = apply_exp_transform
  )
  
  cat("\n--- Function Finished ---\n")
  list(
    settings = list(h = h, pi_level = pi_level, alpha = pi_alpha,
                    exp_transform_applied = apply_exp_transform, sd_method = sd_fts_type),
    optimization_results = list(rolling_only_params = res_roll,
                                expanding_only_params = res_expand,
                                optimal_beta_model = res_optimal),
    test_phase_performance = test_results,
    # Always include RAW test forecasts for convenience
    test_phase_forecasts = c(list(actuals_raw = test_data_raw_actuals,
                                  actuals_smooth = test_data_smooth_actuals),
                             test_forecast_list_raw)
  )
}

# =====================================================================
# 7) EXAMPLE USAGE
# =====================================================================
# Example objects via Female mortality from AUS when h = 1

demog_list$AUS$Female_smooth$rate[[1]] -> smooth_mort_f
demog_list$AUS$Female_raw$rate[[1]]    -> mort_f

final_results <- forecast_with_interval_optimization(
  data_smooth         = log(smooth_mort_f),  # log SMOOTH matrix
  data_raw            = log(mort_f),              # log RAW matrix
  end_train_index     = 61,
  h                   = 1,
  num_components      = 6,
  pi_level            = 0.80,                # desired coverage
  apply_exp_transform = FALSE,                # TRUE => evaluate on RAW; FALSE => SMOOTH/log
  sd_fts_type         = "coordinate",
  seed                = 42,
  objective           = "cpd"                # <-- "cpd" or "interval_score"
)


