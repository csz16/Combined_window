# =====================================================================
# LC_APC_interval_forecast_validated_SD.R
# ---------------------------------------------------------------------
# Validated SD interval forecasts using StMoMo LC/APC point forecasts.
#
# This file keeps the original FTS/SD rolling, expanding, and equal-weight combined structures,
# but replaces the FTS point-forecast engine with StMoMo LC and APC.
#
# Main differences from LC_APC_interval_forecast_3modes.R:
#   * validation set is restored;
#   * SD vectors and tuning parameters are estimated on validation residuals;
#   * test intervals use validation-tuned widths;
#   * equal-weight combined forecasts are included;
#   * h is a user-supplied forecast horizon;
#   * example usage demonstrates h = 1:19 by looping over the single-h
#     wrapper.
#
# Equal-weight combined forecasts use beta = 0.5 and combine the
# validation-tuned rolling/expanding interval half-widths. No optimal-beta
# combined model is used.
#
# Required helper file:
#   LC_APC_point_forecast_3modes.R
# must define:
#   stmomo_fit_and_forecast_all_h()
#   resolve_ages_years_stmomo()
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(StMoMo)
  library(forecast)
  library(ftsa)
})

# =====================================================================
# 0) DEPENDENCY CHECK
# =====================================================================
.ensure_stmomo_point_forecast_helpers <- function(point_forecast_file = "LC_APC_point_forecast_3modes.R") {
  if (!exists("stmomo_fit_and_forecast_all_h", mode = "function") ||
      !exists("resolve_ages_years_stmomo",     mode = "function")) {
    if (file.exists(point_forecast_file)) {
      source(point_forecast_file)
    }
  }

  if (!exists("stmomo_fit_and_forecast_all_h", mode = "function")) {
    stop("Missing stmomo_fit_and_forecast_all_h(). Source LC_APC_point_forecast_3modes.R first, or pass point_forecast_file.")
  }
  if (!exists("resolve_ages_years_stmomo", mode = "function")) {
    stop("Missing resolve_ages_years_stmomo(). Source LC_APC_point_forecast_3modes.R first, or pass point_forecast_file.")
  }

  invisible(TRUE)
}

# =====================================================================
# 1) LC/APC FORECAST GENERATOR
# =====================================================================
# Returns the h-step-ahead forecast vector from a StMoMo LC/APC fit.
# The rolling/expanding schemes are controlled only by start_idx/end_idx.
.stmomo_forecast_iteration <- function(rate_mat_raw,
                                       pop_mat,
                                       start_idx,
                                       end_idx,
                                       h,
                                       model_type,
                                       log_output,
                                       ages_r,
                                       years_r,
                                       cohort_clip = 3L) {
  fc <- tryCatch(
    stmomo_fit_and_forecast_all_h(
      rate_mat_raw = rate_mat_raw,
      pop_mat      = pop_mat,
      start_idx    = start_idx,
      end_idx      = end_idx,
      hmax         = h,
      model_type   = model_type,
      log_output   = log_output,
      ages_num     = ages_r,
      years_num    = years_r,
      cohort_clip  = cohort_clip
    ),
    error = function(e) NULL
  )

  if (is.null(fc)) return(rep(NA_real_, nrow(rate_mat_raw)))
  if (is.vector(fc)) fc <- matrix(fc, ncol = 1L)
  fc <- as.matrix(fc)

  if (nrow(fc) != nrow(rate_mat_raw) && ncol(fc) == nrow(rate_mat_raw)) {
    fc <- t(fc)
  }
  if (nrow(fc) != nrow(rate_mat_raw) || h > ncol(fc)) {
    return(rep(NA_real_, nrow(rate_mat_raw)))
  }

  as.numeric(fc[, h, drop = TRUE])
}

.build_stmomo_forecast_matrix <- function(rate_mat_raw,
                                          pop_mat,
                                          target_indices,
                                          h,
                                          mode = c("rolling", "expanding"),
                                          train_len,
                                          model_type,
                                          log_output,
                                          ages_r,
                                          years_r,
                                          cohort_clip = 3L,
                                          verbose = FALSE) {
  mode <- match.arg(mode)
  n_age <- nrow(rate_mat_raw)
  n_col <- length(target_indices)
  out <- matrix(NA_real_, nrow = n_age, ncol = n_col)

  for (i in seq_len(n_col)) {
    end_idx <- target_indices[i] - h
    start_idx <- if (mode == "rolling") {
      max(1L, end_idx - train_len + 1L)
    } else {
      1L
    }

    out[, i] <- .stmomo_forecast_iteration(
      rate_mat_raw = rate_mat_raw,
      pop_mat      = pop_mat,
      start_idx    = start_idx,
      end_idx      = end_idx,
      h            = h,
      model_type   = model_type,
      log_output   = log_output,
      ages_r       = ages_r,
      years_r      = years_r,
      cohort_clip  = cohort_clip
    )

    if (isTRUE(verbose)) {
      cat(sprintf("    %s/%s forecast %d/%d complete\n", model_type, mode, i, n_col))
    }
  }

  out
}

# =====================================================================
# 2) METRICS
# =====================================================================
# Interval score on the evaluation scale. alpha is the miscoverage level,
# e.g. alpha = 0.20 for an 80% prediction interval.
interval_score <- function(holdout, lb, ub, alpha) {
  holdout <- as.vector(holdout)
  lb      <- as.vector(lb)
  ub      <- as.vector(ub)

  valid <- is.finite(holdout) & is.finite(lb) & is.finite(ub)
  if (!any(valid)) return(c(ECP = NA_real_, CPD = NA_real_, mean_score = Inf))

  holdout <- holdout[valid]
  lb      <- lb[valid]
  ub      <- ub[valid]

  below <- holdout < lb
  above <- holdout > ub

  ecp <- 1 - (sum(below) + sum(above)) / length(holdout)
  cpd <- abs(ecp - (1 - alpha))

  score <- (ub - lb) +
    (2 / alpha) * ((lb - holdout) * below + (holdout - ub) * above)

  c(ECP = ecp, CPD = cpd, mean_score = mean(score))
}

# Build an fts object for ftsa::sd.fts across multiple ftsa versions.
# Some installed versions do not export fts(), so ftsa::fts(...) can fail even
# though sd.fts() is available. This helper tries the exported constructor,
# then the internal constructor, and finally a minimal classed object.
.make_fts_for_sd <- function(y_matrix) {
  x <- 0:(nrow(y_matrix) - 1L)

  if (exists("fts", mode = "function")) {
    out <- tryCatch(fts(x = x, y = y_matrix), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  ns <- asNamespace("ftsa")
  if (exists("fts", envir = ns, mode = "function", inherits = FALSE)) {
    f <- get("fts", envir = ns)
    out <- tryCatch(f(x = x, y = y_matrix), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  out <- list(x = x, y = y_matrix)
  class(out) <- "fts"
  out
}

# Age-specific SD via ftsa::sd.fts, with a row-SD fallback.
calculate_fts_sd <- function(residual_matrix, sd_fts_type = "coordinate") {
  residual_matrix <- as.matrix(residual_matrix)
  if (any(!is.finite(residual_matrix))) residual_matrix[!is.finite(residual_matrix)] <- NA

  resi_fts <- .make_fts_for_sd(residual_matrix)

  sd_vector <- tryCatch(
    {
      ans <- ftsa::sd.fts(resi_fts, method = sd_fts_type)
      if (is.list(ans) && !is.null(ans$y)) ans$y else ans
    },
    error = function(e) {
      warning("ftsa::sd.fts failed: ", conditionMessage(e),
              ". Falling back to row SDs (na.rm = TRUE).")
      apply(residual_matrix, 1L, sd, na.rm = TRUE)
    }
  )

  sd_vector <- as.numeric(sd_vector)
  sd_vector[!is.finite(sd_vector)] <- 0
  sd_vector
}

# Robust console printing for long result tables. forecast_stmomo_interval_3modes()
# can return a base data.frame in some package/version combinations; passing
# n = Inf to base print.data.frame is interpreted incorrectly and can trigger
# "invalid 'na.print' specification". Coerce to tibble first.
.print_long_table <- function(x) {
  if (requireNamespace("tibble", quietly = TRUE)) {
    print(tibble::as_tibble(x), n = Inf, width = Inf)
  } else {
    print(x)
  }
}

.point_metrics <- function(pred, actual) {
  valid <- is.finite(actual) & is.finite(pred)
  if (!any(valid)) return(list(RMSE = NA_real_, MAE = NA_real_))

  y  <- as.vector(actual[valid])
  yh <- as.vector(pred[valid])

  list(
    RMSE = sqrt(mean((y - yh)^2)),
    MAE  = mean(abs(y - yh))
  )
}

# =====================================================================
# 3) TUNING: SD MULTIPLIER FROM VALIDATION RESIDUALS
# =====================================================================
tune_para_objective_fn <- function(par,
                                   forecast_matrix,
                                   validation_data,
                                   sd_vector,
                                   alpha_level) {
  tune_para <- as.numeric(par[1L])
  if (!is.finite(tune_para) || tune_para < 0) return(Inf)

  residual_matrix <- validation_data - forecast_matrix
  n_age <- nrow(residual_matrix)
  n_col <- ncol(residual_matrix)

  sd_vector <- as.numeric(sd_vector)
  if (length(sd_vector) != n_age) sd_vector <- rep_len(sd_vector, n_age)
  sd_vector[!is.finite(sd_vector)] <- 0

  lower_mat <- matrix(-tune_para * sd_vector, nrow = n_age, ncol = n_col)
  upper_mat <- matrix( tune_para * sd_vector, nrow = n_age, ncol = n_col)

  valid <- is.finite(residual_matrix) & is.finite(lower_mat) & is.finite(upper_mat)
  if (!any(valid)) return(Inf)

  covered <- (residual_matrix >= lower_mat) & (residual_matrix <= upper_mat)
  ecp <- sum(covered[valid], na.rm = TRUE) / sum(valid)

  abs(ecp - alpha_level)
}

# Route scale exactly as in the original FTS/SD function:
#   apply_exp_transform = TRUE  -> exp(forecasts), compare to data_raw
#   apply_exp_transform = FALSE -> forecasts unchanged, compare to data_smooth
optimize_tuning_parameter <- function(forecast_matrix_val,
                                      validation_data_raw,
                                      validation_data_smooth,
                                      alpha_level,
                                      apply_exp_transform = TRUE,
                                      sd_fts_type = "coordinate",
                                      seed = 42L) {
  if (isTRUE(apply_exp_transform)) {
    forecasts_for_opt <- exp(forecast_matrix_val)
    actuals_for_opt   <- validation_data_raw
  } else {
    forecasts_for_opt <- forecast_matrix_val
    actuals_for_opt   <- validation_data_smooth
  }

  residuals <- actuals_for_opt - forecasts_for_opt
  sd_vector <- calculate_fts_sd(residuals, sd_fts_type = sd_fts_type)

  set.seed(seed)
  dum_1 <- suppressWarnings(optim(
    par = 2,
    fn = tune_para_objective_fn,
    method = "Nelder-Mead",
    forecast_matrix = forecasts_for_opt,
    validation_data = actuals_for_opt,
    sd_vector = sd_vector,
    alpha_level = alpha_level
  ))

  dum_2 <- optim(
    par = 2,
    fn = tune_para_objective_fn,
    method = "Brent",
    lower = 0,
    upper = 1e3,
    forecast_matrix = forecasts_for_opt,
    validation_data = actuals_for_opt,
    sd_vector = sd_vector,
    alpha_level = alpha_level
  )

  dum_3 <- optim(
    par = 2,
    fn = tune_para_objective_fn,
    method = "L-BFGS-B",
    lower = 0,
    upper = 1e3,
    forecast_matrix = forecasts_for_opt,
    validation_data = actuals_for_opt,
    sd_vector = sd_vector,
    alpha_level = alpha_level
  )

  dum_4 <- optimise(
    f = tune_para_objective_fn,
    interval = c(0, 1e3),
    forecast_matrix = forecasts_for_opt,
    validation_data = actuals_for_opt,
    sd_vector = sd_vector,
    alpha_level = alpha_level
  )

  vals <- c(dum_1$value, dum_2$value, dum_3$value, dum_4$objective)
  k <- which.min(vals)
  best_par <- if (k == 1L) {
    as.numeric(dum_1$par)
  } else if (k == 2L) {
    as.numeric(dum_2$par)
  } else if (k == 3L) {
    as.numeric(dum_3$par)
  } else {
    as.numeric(dum_4$minimum)
  }

  list(
    tune_para       = best_par,
    sd_vector       = sd_vector,
    objective_value = vals[k],
    optimizer_index = k
  )
}

# =====================================================================
# 4) TEST EVALUATION
# =====================================================================
compare_test_performance <- function(forecast_list,
                                     actual_test_data,
                                     pi_alpha,
                                     apply_exp_transform,
                                     verbose = TRUE) {
  method_names <- names(forecast_list)
  results <- data.frame()

  if (isTRUE(verbose)) {
    scale_msg <- if (apply_exp_transform) "RAW (exp/additive)" else "SMOOTH/log or original scale (additive)"
    cat(sprintf("Comparing performance on test set (%s)...\n", scale_msg))
  }

  for (method in method_names) {
    pf <- forecast_list[[method]]

    valid_idx <- is.finite(actual_test_data) & is.finite(pf$mean)
    actuals <- as.vector(actual_test_data[valid_idx])
    preds   <- as.vector(pf$mean[valid_idx])

    rmse <- if (length(actuals) > 0L) sqrt(mean((actuals - preds)^2)) else NA_real_
    mae  <- if (length(actuals) > 0L) mean(abs(actuals - preds)) else NA_real_

    ecp <- NA_real_
    cpd <- NA_real_
    mean_score <- NA_real_

    if (!is.null(pf$lower) && !is.null(pf$upper)) {
      m <- interval_score(
        holdout = actual_test_data,
        lb      = pf$lower,
        ub      = pf$upper,
        alpha   = pi_alpha
      )
      ecp        <- unname(m["ECP"])
      cpd        <- unname(m["CPD"])
      mean_score <- unname(m["mean_score"])
    }

    results <- rbind(
      results,
      data.frame(
        Method = method,
        RMSE = rmse,
        MAE = mae,
        ECP = ecp,
        CPD = cpd,
        Mean_Interval_Score = mean_score,
        Target_Coverage = 1 - pi_alpha,
        Alpha = pi_alpha,
        Exp_Transform = apply_exp_transform,
        row.names = NULL
      )
    )
  }

  results
}

.exp_forecast_list <- function(forecast_list) {
  lapply(forecast_list, function(pf) {
    list(
      mean  = exp(pf$mean),
      lower = exp(pf$lower),
      upper = exp(pf$upper)
    )
  })
}

# =====================================================================
# 5) MAIN SINGLE-H WRAPPER
# =====================================================================
forecast_with_interval_optimization_stmomo <- function(rate_mat_raw,
                                                       pop_mat,
                                                       data_smooth,
                                                       data_raw,
                                                       end_train_index,
                                                       h,
                                                       model_types = c("LC", "APC"),
                                                       log_output = TRUE,
                                                       apply_exp_transform = TRUE,
                                                       pi_level = 0.80,
                                                       sd_fts_type = "coordinate",
                                                       ages = NULL,
                                                       years = NULL,
                                                       cohort_clip = 3L,
                                                       seed = 42L,
                                                       point_forecast_file = "LC_APC_point_forecast_3modes.R",
                                                       verbose = TRUE) {
  .ensure_stmomo_point_forecast_helpers(point_forecast_file = point_forecast_file)

  if (!isTRUE(is.matrix(rate_mat_raw))) rate_mat_raw <- as.matrix(rate_mat_raw)
  if (!isTRUE(is.matrix(pop_mat)))      pop_mat      <- as.matrix(pop_mat)
  if (!isTRUE(is.matrix(data_smooth)))  data_smooth  <- as.matrix(data_smooth)
  if (!isTRUE(is.matrix(data_raw)))     data_raw     <- as.matrix(data_raw)

  if (!all(dim(data_smooth) == dim(data_raw))) {
    stop("data_smooth and data_raw must have the same dimensions.")
  }
  if (nrow(rate_mat_raw) != nrow(data_smooth) || ncol(rate_mat_raw) != ncol(data_smooth)) {
    stop("rate_mat_raw and data_smooth/data_raw must have the same dimensions.")
  }
  if (nrow(pop_mat) != nrow(rate_mat_raw) || ncol(pop_mat) != ncol(rate_mat_raw)) {
    stop("pop_mat and rate_mat_raw must have the same dimensions.")
  }
  if (isTRUE(apply_exp_transform) && !isTRUE(log_output)) {
    stop("apply_exp_transform = TRUE requires log_output = TRUE.")
  }

  h <- as.integer(h)
  if (length(h) != 1L || !is.finite(h) || h < 1L) stop("h must be a positive integer scalar.")

  model_types <- toupper(as.character(model_types))
  bad_models <- setdiff(model_types, c("LC", "APC"))
  if (length(bad_models) > 0L) {
    stop("model_types must contain only 'LC' and/or 'APC'.")
  }

  pi_alpha <- 1 - pi_level
  if (!is.finite(pi_alpha) || pi_alpha <= 0 || pi_alpha >= 1) {
    stop("pi_level must be between 0 and 1.")
  }

  set.seed(seed)

  num_rows    <- nrow(data_smooth)
  num_periods <- ncol(data_smooth)
  train_len   <- end_train_index

  if (end_train_index < 1L || end_train_index >= num_periods) {
    stop("end_train_index must be between 1 and ncol(data_smooth) - 1.")
  }
  if (end_train_index + h > num_periods) {
    stop("Not enough columns after end_train_index for this h.")
  }

  # -------------------------------------------------------------------
  # Validation/test split copied from the original FTS/SD implementation.
  # -------------------------------------------------------------------
  available_all <- seq.int(end_train_index + h, num_periods)
  if (length(available_all) < 2L) stop("Not enough data for validation and test sets.")

  total_start <- num_periods - end_train_index
  base_pairs  <- floor(total_start / 2L)
  pairs_h     <- max(base_pairs - (h - 1L), 1L)
  val_len     <- pairs_h
  test_len    <- pairs_h

  val_pos <- seq_len(val_len)
  test_start_pos <- val_len + h
  test_pos <- seq.int(test_start_pos, test_start_pos + test_len - 1L)

  if (max(test_pos) > length(available_all)) {
    stop(sprintf(
      "Split overflow: need %d positions, but only %d available (h = %d).",
      max(test_pos), length(available_all), h
    ))
  }

  validation_target_indices <- available_all[val_pos]
  test_target_indices       <- available_all[test_pos]

  if (isTRUE(verbose)) {
    cat(sprintf(
      "Split (h=%d): total_start=%d | available=%d | val=%d (%d..%d) | test=%d (%d..%d) | gap=%d\n",
      h, total_start, length(available_all),
      val_len,  min(validation_target_indices), max(validation_target_indices),
      test_len, min(test_target_indices),       max(test_target_indices),
      h - 1L
    ))
  }

  validation_data_raw_actuals    <- data_raw[,    validation_target_indices, drop = FALSE]
  validation_data_smooth_actuals <- data_smooth[, validation_target_indices, drop = FALSE]
  test_data_raw_actuals          <- data_raw[,    test_target_indices,       drop = FALSE]
  test_data_smooth_actuals       <- data_smooth[, test_target_indices,       drop = FALSE]

  ay      <- resolve_ages_years_stmomo(rate_mat_raw, ages = ages, years = years)
  ages_r  <- ay$ages
  years_r <- ay$years

  optimization_by_model <- list()
  performance_by_model  <- list()
  forecasts_eval_by_model <- list()
  forecasts_raw_by_model  <- list()

  for (model_type in model_types) {
    if (isTRUE(verbose)) cat(sprintf("\n--- %s: validation forecasts ---\n", model_type))

    forecast_rates_rolling_val <- .build_stmomo_forecast_matrix(
      rate_mat_raw   = rate_mat_raw,
      pop_mat        = pop_mat,
      target_indices = validation_target_indices,
      h              = h,
      mode           = "rolling",
      train_len      = train_len,
      model_type     = model_type,
      log_output     = log_output,
      ages_r         = ages_r,
      years_r        = years_r,
      cohort_clip    = cohort_clip,
      verbose        = FALSE
    )

    forecast_rates_expanding_val <- .build_stmomo_forecast_matrix(
      rate_mat_raw   = rate_mat_raw,
      pop_mat        = pop_mat,
      target_indices = validation_target_indices,
      h              = h,
      mode           = "expanding",
      train_len      = train_len,
      model_type     = model_type,
      log_output     = log_output,
      ages_r         = ages_r,
      years_r        = years_r,
      cohort_clip    = cohort_clip,
      verbose        = FALSE
    )

    if (isTRUE(verbose)) cat(sprintf("--- %s: optimizing SD tuning parameters ---\n", model_type))

    res_roll <- optimize_tuning_parameter(
      forecast_matrix_val    = forecast_rates_rolling_val,
      validation_data_raw    = validation_data_raw_actuals,
      validation_data_smooth = validation_data_smooth_actuals,
      alpha_level            = pi_level,
      apply_exp_transform    = apply_exp_transform,
      sd_fts_type            = sd_fts_type,
      seed                   = seed
    )

    res_expand <- optimize_tuning_parameter(
      forecast_matrix_val    = forecast_rates_expanding_val,
      validation_data_raw    = validation_data_raw_actuals,
      validation_data_smooth = validation_data_smooth_actuals,
      alpha_level            = pi_level,
      apply_exp_transform    = apply_exp_transform,
      sd_fts_type            = sd_fts_type,
      seed                   = seed
    )

    if (isTRUE(verbose)) {
      cat(sprintf(
        "%s tune_para: rolling = %.4f | expanding = %.4f | combined_equal effective = %.4f\n",
        model_type,
        res_roll$tune_para,
        res_expand$tune_para,
        0.5 * (res_roll$tune_para + res_expand$tune_para)
      ))
      cat(sprintf("--- %s: test forecasts and validation-derived intervals ---\n", model_type))
    }

    fc_rolling_test <- .build_stmomo_forecast_matrix(
      rate_mat_raw   = rate_mat_raw,
      pop_mat        = pop_mat,
      target_indices = test_target_indices,
      h              = h,
      mode           = "rolling",
      train_len      = train_len,
      model_type     = model_type,
      log_output     = log_output,
      ages_r         = ages_r,
      years_r        = years_r,
      cohort_clip    = cohort_clip,
      verbose        = FALSE
    )

    fc_expanding_test <- .build_stmomo_forecast_matrix(
      rate_mat_raw   = rate_mat_raw,
      pop_mat        = pop_mat,
      target_indices = test_target_indices,
      h              = h,
      mode           = "expanding",
      train_len      = train_len,
      model_type     = model_type,
      log_output     = log_output,
      ages_r         = ages_r,
      years_r        = years_r,
      cohort_clip    = cohort_clip,
      verbose        = FALSE
    )

    if (isTRUE(apply_exp_transform)) {
      roll_mean <- exp(fc_rolling_test)
      expd_mean <- exp(fc_expanding_test)

      width_roll <- res_roll$tune_para   * res_roll$sd_vector
      width_expd <- res_expand$tune_para * res_expand$sd_vector

      actuals_for_eval <- test_data_raw_actuals
    } else {
      roll_mean <- fc_rolling_test
      expd_mean <- fc_expanding_test

      width_roll <- res_roll$tune_para   * res_roll$sd_vector
      width_expd <- res_expand$tune_para * res_expand$sd_vector

      actuals_for_eval <- test_data_smooth_actuals
    }

    equal_beta <- 0.5
    combined_mean  <- equal_beta * roll_mean  + (1 - equal_beta) * expd_mean
    width_combined <- equal_beta * width_roll + (1 - equal_beta) * width_expd

    test_forecast_list_for_eval <- list(
      rolling = list(
        mean  = roll_mean,
        lower = sweep(roll_mean, 1L, width_roll, "-"),
        upper = sweep(roll_mean, 1L, width_roll, "+")
      ),
      expanding = list(
        mean  = expd_mean,
        lower = sweep(expd_mean, 1L, width_expd, "-"),
        upper = sweep(expd_mean, 1L, width_expd, "+")
      ),
      combined_equal = list(
        mean  = combined_mean,
        lower = sweep(combined_mean, 1L, width_combined, "-"),
        upper = sweep(combined_mean, 1L, width_combined, "+")
      )
    )

    perf <- compare_test_performance(
      forecast_list       = test_forecast_list_for_eval,
      actual_test_data    = actuals_for_eval,
      pi_alpha            = pi_alpha,
      apply_exp_transform = apply_exp_transform,
      verbose             = verbose
    )

    perf <- perf %>%
      dplyr::mutate(
        Model = model_type,
        h = h,
        Tune_Para = dplyr::case_when(
          Method == "rolling"        ~ res_roll$tune_para,
          Method == "expanding"      ~ res_expand$tune_para,
          Method == "combined_equal" ~ equal_beta * res_roll$tune_para +
                                        (1 - equal_beta) * res_expand$tune_para,
          TRUE ~ NA_real_
        ),
        SD_Method = sd_fts_type
      ) %>%
      dplyr::select(Model, h, Method, RMSE, MAE, ECP, CPD,
                    Mean_Interval_Score, Target_Coverage, Alpha,
                    Exp_Transform, Tune_Para, SD_Method)

    optimization_by_model[[model_type]] <- list(
      rolling_only_params   = res_roll,
      expanding_only_params = res_expand,
      combined_equal_params = list(
        beta = equal_beta,
        tune_para = equal_beta * res_roll$tune_para +
          (1 - equal_beta) * res_expand$tune_para,
        rolling_tune_para = res_roll$tune_para,
        expanding_tune_para = res_expand$tune_para,
        width_vector = width_combined,
        width_vector_source = "0.5 * rolling_width + 0.5 * expanding_width"
      )
    )

    performance_by_model[[model_type]] <- perf
    forecasts_eval_by_model[[model_type]] <- test_forecast_list_for_eval

    forecasts_raw_by_model[[model_type]] <- if (isTRUE(apply_exp_transform)) {
      test_forecast_list_for_eval
    } else if (isTRUE(log_output)) {
      .exp_forecast_list(test_forecast_list_for_eval)
    } else {
      test_forecast_list_for_eval
    }
  }

  test_results <- dplyr::bind_rows(performance_by_model)

  if (isTRUE(verbose)) cat("\n--- Function Finished ---\n")

  list(
    settings = list(
      h = h,
      model_types = model_types,
      pi_level = pi_level,
      alpha = pi_alpha,
      exp_transform_applied = apply_exp_transform,
      log_output = log_output,
      sd_method = sd_fts_type,
      end_train_index = end_train_index,
      validation_target_indices = validation_target_indices,
      test_target_indices = test_target_indices,
      rolling_window_length = train_len
    ),
    optimization_results = optimization_by_model,
    test_phase_performance = test_results,
    test_phase_forecasts = list(
      actuals_raw = test_data_raw_actuals,
      actuals_smooth = test_data_smooth_actuals,
      evaluation_scale = forecasts_eval_by_model,
      raw_scale = forecasts_raw_by_model
    )
  )
}

# Backwards-compatible alias with a shorter StMoMo-specific name.
forecast_stmomo_interval_validated_SD <- forecast_with_interval_optimization_stmomo

# =====================================================================
# 6) MULTI-HORIZON HELPER
# =====================================================================
forecast_with_interval_optimization_stmomo_horizons <- function(rate_mat_raw,
                                                                pop_mat,
                                                                data_smooth,
                                                                data_raw,
                                                                end_train_index,
                                                                h_values = 1:19,
                                                                ...) {
  out <- lapply(as.integer(h_values), function(h_i) {
    forecast_with_interval_optimization_stmomo(
      rate_mat_raw     = rate_mat_raw,
      pop_mat          = pop_mat,
      data_smooth      = data_smooth,
      data_raw         = data_raw,
      end_train_index  = end_train_index,
      h                = h_i,
      ...
    )
  })
  names(out) <- paste0("h", as.integer(h_values))
  out
}

stmomo_interval_performance_table <- function(results_by_h) {
  dplyr::bind_rows(lapply(results_by_h, function(x) x$test_phase_performance))
}

# Compatibility helper that returns the long table shape used by the
# previous LC/APC 3-mode interval file. Rows include rolling, expanding,
# and combined_equal. The combined_equal row uses beta = 0.5 only.
forecast_stmomo_interval_3modes <- function(rate_mat_raw,
                                            pop_mat,
                                            data_smooth,
                                            data_raw,
                                            end_train_index,
                                            h_values = 1:19,
                                            model_types = c("LC", "APC"),
                                            log_output = TRUE,
                                            apply_exp_transform = TRUE,
                                            pi_level = 0.80,
                                            sd_fts_type = "coordinate",
                                            ages = NULL,
                                            years = NULL,
                                            cohort_clip = 3L,
                                            seed = 42L,
                                            point_forecast_file = "LC_APC_point_forecast_3modes.R",
                                            verbose = TRUE) {
  results_by_h <- forecast_with_interval_optimization_stmomo_horizons(
    rate_mat_raw        = rate_mat_raw,
    pop_mat             = pop_mat,
    data_smooth         = data_smooth,
    data_raw            = data_raw,
    end_train_index     = end_train_index,
    h_values            = h_values,
    model_types         = model_types,
    log_output          = log_output,
    apply_exp_transform = apply_exp_transform,
    pi_level            = pi_level,
    sd_fts_type         = sd_fts_type,
    ages                = ages,
    years               = years,
    cohort_clip         = cohort_clip,
    seed                = seed,
    point_forecast_file = point_forecast_file,
    verbose             = verbose
  )

  stmomo_interval_performance_table(results_by_h) %>%
    dplyr::rename(method = Model, mode = Method, Mean_IS = Mean_Interval_Score) %>%
    dplyr::mutate(status = "ok", error = NA_character_, pi_level = Target_Coverage) %>%
    dplyr::select(h, method, mode, RMSE, MAE, ECP, CPD, Mean_IS,
                  pi_level, Alpha, Exp_Transform, Tune_Para, SD_Method,
                  status, error)
}

# =====================================================================
# 7) SPLIT INSPECTOR
# =====================================================================
inspect_split_stmomo_interval_validated_SD <- function(num_periods,
                                                       end_train_index,
                                                       h_values = 1:19) {
  total_start <- num_periods - end_train_index
  base_pairs  <- as.integer(floor(total_start / 2L))

  cat(sprintf("num_periods     = %d\n", num_periods))
  cat(sprintf("end_train_index = %d\n", end_train_index))
  cat(sprintf("total_start     = %d\n", total_start))
  cat(sprintf("base_pairs      = %d\n\n", base_pairs))

  cat(sprintf("%-4s %-6s %-18s %-18s %-5s\n", "h", "n_each", "validation targets", "test targets", "gap"))
  cat(sprintf("%-4s %-6s %-18s %-18s %-5s\n", "---", "------", "------------------", "------------", "---"))

  for (h in as.integer(h_values)) {
    if (end_train_index + h > num_periods) {
      cat(sprintf("%-4d %-6s %-18s %-18s %-5s\n", h, "0", "(none)", "(none)", "--"))
      next
    }

    available_all <- seq.int(end_train_index + h, num_periods)
    pairs_h <- max(base_pairs - (h - 1L), 1L)
    val_pos <- seq_len(pairs_h)
    test_start_pos <- pairs_h + h
    test_pos <- seq.int(test_start_pos, test_start_pos + pairs_h - 1L)

    if (max(test_pos) > length(available_all)) {
      cat(sprintf("%-4d %-6s %-18s %-18s %-5s\n", h, "--", "overflow", "overflow", "--"))
      next
    }

    val_idx  <- available_all[val_pos]
    test_idx <- available_all[test_pos]

    cat(sprintf(
      "%-4d %-6d %3d..%-13d %3d..%-13d %-5d\n",
      h, pairs_h, min(val_idx), max(val_idx), min(test_idx), max(test_idx), h - 1L
    ))
  }

  invisible(NULL)
}

# =====================================================================
# 8) EXAMPLE USAGE
# =====================================================================
# Example 1: single horizon, output format similar to the original FTS/SD
# implementation.
#
# demog_list$AUS$Female_smooth$rate[[1]] -> smooth_mort_f
# demog_list$AUS$Female_raw$rate[[1]]    -> mort_f
# demog_list$AUS$Female_raw$pop[[1]]     -> pop_mort_f
#
# final_results_h1 <- forecast_with_interval_optimization_stmomo(
#   rate_mat_raw        = mort_f,
#   pop_mat             = pop_mort_f,
#   data_smooth         = log(smooth_mort_f),
#   data_raw            = log(mort_f),
#   end_train_index     = 61,
#   h                   = 1,
#   model_types         = c("LC", "APC"),
#   log_output          = TRUE,
#   apply_exp_transform = FALSE,
#   pi_level            = 0.80,
#   sd_fts_type         = "coordinate",
#   cohort_clip         = 3L,
#   seed                = 42L
# )
#
# final_results_h1$test_phase_performance
# final_results_h1$optimization_results$LC$rolling_only_params$tune_para
# final_results_h1$optimization_results$LC$expanding_only_params$tune_para
# final_results_h1$optimization_results$LC$combined_equal_params$tune_para
# final_results_h1$optimization_results$APC$combined_equal_params$tune_para
#
# Example 2: apply the same single-h wrapper over h = 1:19.
#
# results_h_1_to_19 <- forecast_with_interval_optimization_stmomo_horizons(
#   rate_mat_raw        = mort_f,
#   pop_mat             = pop_mort_f,
#   data_smooth         = log(smooth_mort_f),
#   data_raw            = log(mort_f),
#   end_train_index     = 61,
#   h_values            = 1:19,
#   model_types         = c("LC", "APC"),
#   log_output          = TRUE,
#   apply_exp_transform = FALSE,
#   pi_level            = 0.80,
#   sd_fts_type         = "coordinate",
#   cohort_clip         = 3L,
#   seed                = 42L,
#   verbose             = FALSE
# )
#
# performance_h_1_to_19 <- stmomo_interval_performance_table(results_h_1_to_19)
# performance_h_1_to_19
#
# Optional compatibility output matching the previous long-table LC/APC file:
#
# long_table_h_1_to_19 <- forecast_stmomo_interval_3modes(
#   rate_mat_raw        = mort_f,
#   pop_mat             = pop_mort_f,
#   data_smooth         = log(smooth_mort_f),
#   data_raw            = log(mort_f),
#   end_train_index     = 61,
#   h_values            = 1:19,
#   model_types         = c("LC", "APC"),
#   log_output          = TRUE,
#   apply_exp_transform = FALSE,
#   pi_level            = 0.80,
#   sd_fts_type         = "coordinate",
#   cohort_clip         = 3L,
#   seed                = 42L,
#   verbose             = FALSE
# )
# =====================================================================

# =====================================================================
# 9) OPTIONAL EXAMPLE AND FULL-PIPELINE DRIVERS
# =====================================================================
# These drivers retain the calling style of the previous LC/APC interval
# file, but now use the validated SD workflow with three schemes:
# rolling, expanding, and combined_equal.

run_example_stmomo_interval_validated_SD <- function(pi_level = 0.80,
                                                     h_values = 1:19,
                                                     model_types = c("LC", "APC")) {
  cat("\n========== EXAMPLE: VALIDATED SD LC/APC INTERVALS ==========" , "\n")

  # ---- Mortality: AUS / Female ----
  if (!file.exists("demog_list.rds")) {
    cat("[skip] demog_list.rds not found.\n")
  } else {
    demog_list <- readRDS("demog_list.rds")
    if (is.null(demog_list$AUS)) {
      cat("[skip] demog_list$AUS not found.\n")
    } else {
      e <- demog_list$AUS
      rate_raw    <- e$Female_raw$rate[[1]]
      rate_smooth <- e$Female_smooth$rate[[1]]
      pop         <- e$Female_raw$pop[[1]]
      nT <- ncol(rate_raw)
      end_train_index <- nT - 40L

      cat(sprintf("\n--- Mortality: AUS / Female (nT=%d, end_train=%d) ---\n",
                  nT, end_train_index))

      mort_df <- forecast_stmomo_interval_3modes(
        rate_mat_raw        = rate_raw,
        pop_mat             = pop,
        data_smooth         = log(rate_smooth),
        data_raw            = log(rate_raw),
        end_train_index     = end_train_index,
        h_values            = h_values,
        model_types         = model_types,
        log_output          = TRUE,
        apply_exp_transform = FALSE,
        pi_level            = pi_level,
        cohort_clip         = 3L,
        seed                = 42L,
        verbose             = FALSE
      )

      cat("\nSummary by (method, mode), mean over h:\n")
      print(
        mort_df %>% dplyr::filter(status == "ok") %>%
          dplyr::group_by(method, mode) %>%
          dplyr::summarise(RMSE_mean    = mean(RMSE,    na.rm = TRUE),
                           MAE_mean     = mean(MAE,     na.rm = TRUE),
                           ECP_mean     = mean(ECP,     na.rm = TRUE),
                           CPD_mean     = mean(CPD,     na.rm = TRUE),
                           Mean_IS_mean = mean(Mean_IS, na.rm = TRUE),
                           .groups = "drop")
      )

      cat("\nFull long-format table:\n")
      .print_long_table(mort_df)
    }
  }

  # ---- ASFR: CAN ----
  if (!file.exists("demog_list_asfr.rds")) {
    cat("\n[skip] demog_list_asfr.rds not found.\n")
  } else {
    demog_list_asfr <- readRDS("demog_list_asfr.rds")
    if (is.null(demog_list_asfr$CAN)) {
      cat("\n[skip] demog_list_asfr$CAN not found.\n")
    } else {
      e <- demog_list_asfr$CAN
      rate_raw    <- as.matrix(e$raw$rate[[1]])
      rate_smooth <- as.matrix(e$smooth$rate[[1]])
      pop         <- as.matrix(e$raw$pop[[1]])
      nT <- ncol(rate_raw)
      end_train_index <- nT - 40L

      cat(sprintf("\n--- ASFR: CAN (nT=%d, end_train=%d) ---\n",
                  nT, end_train_index))

      asfr_df <- forecast_stmomo_interval_3modes(
        rate_mat_raw        = rate_raw,
        pop_mat             = pop,
        data_smooth         = rate_smooth,
        data_raw            = rate_raw,
        end_train_index     = end_train_index,
        h_values            = h_values,
        model_types         = model_types,
        log_output          = FALSE,
        apply_exp_transform = FALSE,
        pi_level            = pi_level,
        cohort_clip         = 3L,
        seed                = 42L,
        verbose             = FALSE
      )

      cat("\nSummary by (method, mode), mean over h:\n")
      print(
        asfr_df %>% dplyr::filter(status == "ok") %>%
          dplyr::group_by(method, mode) %>%
          dplyr::summarise(RMSE_mean    = mean(RMSE,    na.rm = TRUE),
                           MAE_mean     = mean(MAE,     na.rm = TRUE),
                           ECP_mean     = mean(ECP,     na.rm = TRUE),
                           CPD_mean     = mean(CPD,     na.rm = TRUE),
                           Mean_IS_mean = mean(Mean_IS, na.rm = TRUE),
                           .groups = "drop")
      )

      cat("\nFull long-format table:\n")
      .print_long_table(asfr_df)
    }
  }

  cat("\nTo run the full pipeline: run_all_stmomo_interval_validated_SD()\n")
  invisible(NULL)
}

# Backwards-compatible driver name. It now runs rolling, expanding, and
# equal-weight combined rows.
run_example_3modes_interval <- run_example_stmomo_interval_validated_SD

run_all_stmomo_interval_validated_SD <- function(
  h_values    = 1:19,
  model_types = c("LC", "APC"),
  pi_level    = 0.80,
  sd_fts_type = "coordinate",
  cohort_clip = 3L,
  seed        = 42L,
  out_dir     = "results/interval_validated_SD",
  mort_file   = "results_interval_validated_SD_mortality.rds",
  asfr_file   = "results_interval_validated_SD_asfr.rds",
  all_file    = "results_interval_validated_SD_ALL.rds"
) {
  expected_rows_per_unit <- length(model_types) * 3L * length(h_values)

  is_unit_complete <- function(df, dom, cc, sx) {
    if (!nrow(df)) return(FALSE)
    n <- if (is.na(sx)) {
      sum(df$domain == dom & df$country == cc & is.na(df$sex))
    } else {
      sum(df$domain == dom & df$country == cc & !is.na(df$sex) & df$sex == sx)
    }
    n >= expected_rows_per_unit
  }

  empty_tbl <- tibble::tibble(
    domain   = character(), country = character(), sex   = character(),
    h        = integer(),   method  = character(), mode  = character(),
    RMSE     = numeric(),   MAE     = numeric(),
    ECP      = numeric(),   CPD     = numeric(),   Mean_IS = numeric(),
    pi_level = numeric(),   Alpha   = numeric(),   Exp_Transform = logical(),
    Tune_Para = numeric(),  SD_Method = character(),
    status   = character(), error   = character()
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mort_path <- file.path(out_dir, mort_file)
  asfr_path <- file.path(out_dir, asfr_file)
  all_path  <- file.path(out_dir, all_file)

  results_mort <- if (file.exists(mort_path)) readRDS(mort_path) else empty_tbl
  results_asfr <- if (file.exists(asfr_path)) readRDS(asfr_path) else empty_tbl

  units_all <- list()

  demog_list <- NULL
  if (file.exists("demog_list.rds")) {
    demog_list <- readRDS("demog_list.rds")
    mort_countries <- names(demog_list)[!vapply(demog_list, is.null, logical(1L))]
    for (cc in mort_countries) {
      for (sex_key in c("Female", "Male")) {
        units_all[[length(units_all) + 1L]] <-
          list(domain = "mortality", country = cc, sex = sex_key)
      }
    }
  }

  demog_list_asfr <- NULL
  if (file.exists("demog_list_asfr.rds")) {
    demog_list_asfr <- readRDS("demog_list_asfr.rds")
    asfr_countries <- names(demog_list_asfr)[!vapply(demog_list_asfr, is.null, logical(1L))]
    for (cc in asfr_countries) {
      units_all[[length(units_all) + 1L]] <-
        list(domain = "asfr", country = cc, sex = NA_character_)
    }
  }

  total_units <- length(units_all)
  if (!total_units) {
    cat("No data found: demog_list.rds and demog_list_asfr.rds are missing.\n")
    return(invisible(NULL))
  }

  units_todo <- list()
  units_done_carry <- 0L
  for (u in units_all) {
    df_check <- if (u$domain == "mortality") results_mort else results_asfr
    if (is_unit_complete(df_check, u$domain, u$country, u$sex)) {
      units_done_carry <- units_done_carry + 1L
    } else {
      units_todo[[length(units_todo) + 1L]] <- u
    }
  }

  cat("\n========== FULL RUN: VALIDATED SD LC/APC INTERVALS ==========" , "\n")
  cat(sprintf("pi_level                    : %.2f\n", pi_level))
  cat(sprintf("h_values                    : %s\n", paste(h_values, collapse = ",")))
  cat(sprintf("models                      : %s\n", paste(model_types, collapse = ",")))
  cat("schemes                     : rolling, expanding, combined_equal\n")
  cat(sprintf("rows per unit               : %d\n", expected_rows_per_unit))
  cat(sprintf("total units                 : %d\n", total_units))
  cat(sprintf("already complete            : %d\n", units_done_carry))
  cat(sprintf("to run now                  : %d\n", length(units_todo)))

  if (!length(units_todo)) {
    results_all <- dplyr::bind_rows(results_mort, results_asfr)
    saveRDS(results_all, all_path)
    return(invisible(list(mort = results_mort,
                          asfr = results_asfr,
                          all  = results_all)))
  }

  for (idx in seq_along(units_todo)) {
    u <- units_todo[[idx]]
    sxL <- if (is.na(u$sex)) "" else paste0("/", u$sex)
    cat(sprintf("\n[%d/%d] %s %s%s\n", idx, length(units_todo), u$domain, u$country, sxL))

    if (u$domain == "mortality") {
      entry <- demog_list[[u$country]]
      rk <- paste0(u$sex, "_raw")
      sk <- paste0(u$sex, "_smooth")

      if (is.null(entry[[rk]]) || is.null(entry[[sk]]) ||
          is.null(entry[[rk]]$rate[[1]]) || is.null(entry[[sk]]$rate[[1]]) ||
          is.null(entry[[rk]]$pop[[1]])) {
        cat("  [skip] missing mortality rate/pop data\n")
        next
      }

      rate_raw    <- entry[[rk]]$rate[[1]]
      rate_smooth <- entry[[sk]]$rate[[1]]
      pop         <- entry[[rk]]$pop[[1]]
      nT <- ncol(rate_raw)
      end_train_index <- nT - 40L
      if (end_train_index < 1L) {
        cat("  [skip] nT too small\n")
        next
      }

      df <- forecast_stmomo_interval_3modes(
        rate_mat_raw        = rate_raw,
        pop_mat             = pop,
        data_smooth         = log(rate_smooth),
        data_raw            = log(rate_raw),
        end_train_index     = end_train_index,
        h_values            = h_values,
        model_types         = model_types,
        log_output          = TRUE,
        apply_exp_transform = FALSE,
        pi_level            = pi_level,
        sd_fts_type         = sd_fts_type,
        cohort_clip         = cohort_clip,
        seed                = seed,
        verbose             = FALSE
      ) %>%
        dplyr::mutate(domain = "mortality", country = u$country, sex = u$sex) %>%
        dplyr::select(domain, country, sex, h, method, mode,
                      RMSE, MAE, ECP, CPD, Mean_IS,
                      pi_level, Alpha, Exp_Transform, Tune_Para, SD_Method,
                      status, error)

      # Replace any earlier cached rows for this unit. This prevents old
      # rolling/expanding-only runs from being mixed with the new
      # rolling/expanding/combined_equal output.
      results_mort <- results_mort %>%
        dplyr::filter(!(domain == "mortality" & country == u$country &
                          !is.na(sex) & sex == u$sex))
      results_mort <- dplyr::bind_rows(results_mort, df)
      saveRDS(results_mort, mort_path)
    } else {
      entry <- demog_list_asfr[[u$country]]
      if (is.null(entry$raw$rate[[1]]) || is.null(entry$smooth$rate[[1]]) ||
          is.null(entry$raw$pop[[1]])) {
        cat("  [skip] missing ASFR rate/pop data\n")
        next
      }

      rate_raw    <- as.matrix(entry$raw$rate[[1]])
      rate_smooth <- as.matrix(entry$smooth$rate[[1]])
      pop         <- as.matrix(entry$raw$pop[[1]])
      nT <- ncol(rate_raw)
      end_train_index <- nT - 40L
      if (end_train_index < 1L) {
        cat("  [skip] nT too small\n")
        next
      }

      df <- forecast_stmomo_interval_3modes(
        rate_mat_raw        = rate_raw,
        pop_mat             = pop,
        data_smooth         = rate_smooth,
        data_raw            = rate_raw,
        end_train_index     = end_train_index,
        h_values            = h_values,
        model_types         = model_types,
        log_output          = FALSE,
        apply_exp_transform = FALSE,
        pi_level            = pi_level,
        sd_fts_type         = sd_fts_type,
        cohort_clip         = cohort_clip,
        seed                = seed,
        verbose             = FALSE
      ) %>%
        dplyr::mutate(domain = "asfr", country = u$country, sex = NA_character_) %>%
        dplyr::select(domain, country, sex, h, method, mode,
                      RMSE, MAE, ECP, CPD, Mean_IS,
                      pi_level, Alpha, Exp_Transform, Tune_Para, SD_Method,
                      status, error)

      # Replace any earlier cached rows for this unit. This prevents old
      # rolling/expanding-only runs from being mixed with the new
      # rolling/expanding/combined_equal output.
      results_asfr <- results_asfr %>%
        dplyr::filter(!(domain == "asfr" & country == u$country))
      results_asfr <- dplyr::bind_rows(results_asfr, df)
      saveRDS(results_asfr, asfr_path)
    }
  }

  results_all <- dplyr::bind_rows(results_mort, results_asfr)
  saveRDS(results_all, all_path)

  cat("\n========== SUMMARY ==========" , "\n")
  cat(sprintf("Mortality rows : %d\n", nrow(results_mort)))
  cat(sprintf("ASFR rows      : %d\n", nrow(results_asfr)))
  cat(sprintf("Total rows     : %d\n", nrow(results_all)))
  cat(sprintf("Files saved at : %s\n", out_dir))

  print(
    results_all %>%
      dplyr::filter(status == "ok") %>%
      dplyr::group_by(domain, method, mode) %>%
      dplyr::summarise(n_rows         = dplyr::n(),
                       Mean_IS_median = median(Mean_IS, na.rm = TRUE),
                       ECP_median     = median(ECP,     na.rm = TRUE),
                       CPD_median     = median(CPD,     na.rm = TRUE),
                       .groups = "drop")
  )

  invisible(list(mort = results_mort,
                 asfr = results_asfr,
                 all  = results_all))
}

# Backwards-compatible full-run driver name. It now runs rolling, expanding,
# and equal-weight combined rows.
run_all_3modes_interval <- run_all_stmomo_interval_validated_SD

# =====================================================================
# How to use:
  source("LC_APC_interval_forecast_validated_SD.R")
  run_example_stmomo_interval_validated_SD()
#   run_all_stmomo_interval_validated_SD()
# =====================================================================
