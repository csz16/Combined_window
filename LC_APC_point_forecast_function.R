# =====================================================================
# LC_APC_point_forecast_function.R
# -------------------------------------------------------------------
# Point forecasts for:
#   (1) Lee & Carter (1992):
#          log m_{x,t} = alpha_x + beta_x^{(1)} * kappa_t^{(1)}
#          identifiability:  sum_x beta_x = 1,  sum_t kappa_t = 0.
#   (2) Classical Age-Period-Cohort (APC) model in the StMoMo
#       parameterisation (Table 1 of Villegas, Kaishev & Millossovich
#       2018; as implemented by StMoMo::apc). This is the same APC
#       predictor used e.g. in Currie (2006):
#          log m_{x,t} = alpha_x + kappa_t^{(1)} + gamma_{t-x}^{(0)}
#          identifiability:  sum_t kappa_t = 0,
#                            sum_c gamma_c = 0,
#                            sum_c c * gamma_c = 0.
#       Note: this is NOT the Renshaw-Haberman (2006) cohort extension
#       of Lee-Carter (StMoMo::rh), which has an additional age-modulation
#       beta_x^{(0)} on the cohort term.
#
# Estimation: Poisson(D_{x,t}) with log link and offset log(E_{x,t}),
# as implemented in StMoMo::fit(). For data given as rates (mortality
# m_{x,t} or fertility ASFR f_{x,t}) we reconstruct an integer count
# D_{x,t} = round(rate_{x,t} * E_{x,t}) and use central exposures E_{x,t}.
#
# For APC we additionally multiply the validity weight matrix wxt by
# StMoMo::genWeightMat(ages, years, clip = cohort_clip) to down-weight
# the weakly-identified boundary cohorts (default clip = 3, as in the
# StMoMo vignette). This is critical for stable APC fits and has no
# effect for LC (but we still apply it only to the APC branch).
#
# The train / validation / test split exactly mirrors the FTS-based
# Point_forecast_function.R so the comparison is head-to-head:
#     total_start = ncol(data) - end_train_index
#     base_pairs  = floor(total_start / 2)
#     pairs_h     = max(base_pairs - (h - 1), 1)  =>  val_len = test_len = pairs_h
# =====================================================================

suppressPackageStartupMessages({
  library(StMoMo)
  library(forecast)
  library(dplyr)
})

# Session-level flag set: each key fires a warning at most once per
# R session. Users who want the full stream of warnings back can do
# `stmomo_reset_warn_flags()` before running; users who want total
# silence (having already accepted that sequential-index fallback is
# mathematically safe for LC and APC) can do
# `stmomo_suppress_dimname_warnings()` once at the top of a script.
.stmomo_warn_env <- new.env(parent = emptyenv())
stmomo_reset_warn_flags <- function() {
  rm(list = ls(envir = .stmomo_warn_env, all.names = TRUE),
     envir = .stmomo_warn_env)
  invisible(NULL)
}
stmomo_suppress_dimname_warnings <- function() {
  assign("ages_fallback",  TRUE, envir = .stmomo_warn_env)
  assign("years_fallback", TRUE, envir = .stmomo_warn_env)
  invisible(NULL)
}
.stmomo_warn_once <- function(key, msg) {
  if (!isTRUE(.stmomo_warn_env[[key]])) {
    warning(msg, call. = FALSE)
    assign(key, TRUE, envir = .stmomo_warn_env)
  }
  invisible(NULL)
}

# =====================================================================
# 0) INTERNAL HELPER - resolve ages / years from dimnames with fallback
# =====================================================================
# Resolution order for each of `ages` and `years`:
#   1. Use the explicit argument if provided.
#   2. Otherwise, parse the corresponding dimnames of rate_mat_raw.
#      Labels are cleaned before parsing: leading/trailing whitespace
#      is trimmed, and a trailing "+" (HMD/HFD convention marking an
#      open-ended age group, e.g. "95+", "49+") is stripped.
#      "95+" -> 95, "49+" -> 49. This is label-only; see note below.
#   3. If, after cleaning, labels are still not all parseable, fall
#      back to sequential indexing (0, 1, ..., n-1) and emit a warning
#      at most ONCE PER SESSION per warning key.
#
# Mathematical note: both LC and APC are equivariant under constant
# shifts in age and year labels. For LC, a shift is absorbed into
# alpha_x and beta_x. For APC, cohort equivalence classes depend only
# on c = t - x, so any uniform shift of ages leaves the partition --
# and hence the estimates of gamma_c -- unchanged. ARIMA forecasts of
# kappa_t and gamma_c depend only on the ordering of time/cohort
# indices, not the origin. Consequently:
#   (i)  mapping "95+" -> 95 does not perturb fit or forecast, because
#        the open-ended group becomes simply "the last age row" with
#        a numeric label that does not collide with any other;
#   (ii) the sequential-index fallback, on the rare paths where it
#        triggers, preserves in-sample fit and h-step-ahead forecasts
#        exactly; only the cohort LABELS become nominal rather than
#        calendar-accurate.
.parse_demog_labels <- function(x) {
  if (is.null(x) || length(x) == 0L) return(numeric(0))
  cleaned <- trimws(as.character(x))
  cleaned <- sub("\\+$", "", cleaned)   # open-age marker (HMD/HFD)
  suppressWarnings(as.numeric(cleaned))
}

resolve_ages_years_stmomo <- function(rate_mat_raw, ages = NULL, years = NULL) {
  # --- ages ---
  if (!is.null(ages)) {
    # Accept a numeric or character vector; clean "95+" style markers.
    # If the caller passed a pre-parsed numeric vector that already
    # contains NAs (e.g. `as.numeric(c("0","1",...,"95+"))`), try to
    # recover from the underlying dimnames before giving up.
    ages_num <- .parse_demog_labels(ages)
    if (length(ages_num) != nrow(rate_mat_raw) || anyNA(ages_num)) {
      rn_parsed <- .parse_demog_labels(rownames(rate_mat_raw))
      if (length(rn_parsed) == nrow(rate_mat_raw) && !anyNA(rn_parsed)) {
        ages_num <- rn_parsed
      } else {
        stop("`ages` must be coercible to a numeric vector of length nrow(rate_mat_raw).")
      }
    }
  } else {
    ages_num <- .parse_demog_labels(rownames(rate_mat_raw))
    if (length(ages_num) == 0L || anyNA(ages_num)) {
      .stmomo_warn_once(
        "ages_fallback",
        paste0("Row names of rate_mat_raw are not parseable as numeric ages ",
               "(even after trimming whitespace and stripping a trailing '+'); ",
               "falling back to 0, 1, ..., nrow-1. LC and APC fits are invariant ",
               "under constant shifts of age labels, so numerical results are ",
               "unaffected; cohort labels, however, will be nominal. ",
               "(This warning is emitted once per R session; ",
               "call stmomo_reset_warn_flags() to re-enable it, or ",
               "stmomo_suppress_dimname_warnings() to silence it.)")
      )
      ages_num <- seq_len(nrow(rate_mat_raw)) - 1L
    }
  }
  
  # --- years ---
  if (!is.null(years)) {
    yrs_num <- .parse_demog_labels(years)
    if (length(yrs_num) != ncol(rate_mat_raw) || anyNA(yrs_num)) {
      cn_parsed <- .parse_demog_labels(colnames(rate_mat_raw))
      if (length(cn_parsed) == ncol(rate_mat_raw) && !anyNA(cn_parsed)) {
        yrs_num <- cn_parsed
      } else {
        stop("`years` must be coercible to a numeric vector of length ncol(rate_mat_raw).")
      }
    }
  } else {
    yrs_num <- .parse_demog_labels(colnames(rate_mat_raw))
    if (length(yrs_num) == 0L || anyNA(yrs_num)) {
      .stmomo_warn_once(
        "years_fallback",
        paste0("Column names of rate_mat_raw are not parseable as numeric years ",
               "(even after trimming whitespace and stripping a trailing '+'); ",
               "falling back to 0, 1, ..., ncol-1. LC and APC fits and h-step ",
               "forecasts are invariant under constant shifts of year labels. ",
               "(This warning is emitted once per R session; ",
               "call stmomo_reset_warn_flags() to re-enable it, or ",
               "stmomo_suppress_dimname_warnings() to silence it.)")
      )
      yrs_num <- seq_len(ncol(rate_mat_raw)) - 1L
    }
  }
  
  list(ages = ages_num, years = yrs_num)
}

# =====================================================================
# 1) SINGLE FORECAST ITERATION (one training origin, one horizon h)
# =====================================================================
# Returns a vector of length nrow(rate_mat_raw) with the h-step-ahead
# rate forecast (raw rate scale if log_output = FALSE, log-rate scale
# if log_output = TRUE). On any fit/forecast failure, returns NA's.
#
# Arguments:
#   rate_mat_raw : raw rate matrix m_{x,t} (ages in rows, years in columns).
#   pop_mat      : central exposures E_{x,t} (same shape as rate_mat_raw).
#   start_idx, end_idx, h : training window [start_idx, end_idx], horizon h.
#   model_type   : "LC" or "APC".
#   log_output   : TRUE  -> return log m_{x,t+h}; FALSE -> return m_{x,t+h}.
#   ages, years  : OPTIONAL numeric vectors of ages (length = nrow) and
#                  years (length = ncol). If NULL, they are parsed from
#                  the matrix dimnames and an error is thrown when the
#                  dimnames are not numeric. No zero-based fallback is
#                  used because for APC the cohort index is c = t - x,
#                  so an incorrect age origin shifts every cohort and
#                  yields a materially different fit.
#   cohort_clip  : integer (default 3) passed to StMoMo::genWeightMat().
#                  Only has effect for APC.
stmomo_forecast_iter <- function(rate_mat_raw, pop_mat,
                                 start_idx, end_idx, h,
                                 model_type  = c("LC", "APC"),
                                 log_output  = TRUE,
                                 ages        = NULL,
                                 years       = NULL,
                                 cohort_clip = 3L,
                                 verbose_fit = FALSE) {
  model_type <- match.arg(model_type)
  
  rate_w <- rate_mat_raw[, start_idx:end_idx, drop = FALSE]
  pop_w  <- pop_mat[,     start_idx:end_idx, drop = FALSE]
  
  # --- ages ---
  if (is.null(ages)) {
    ages_num <- suppressWarnings(as.numeric(rownames(rate_mat_raw)))
    if (anyNA(ages_num) || length(ages_num) == 0L) {
      stop("stmomo_forecast_iter(): numeric `ages` vector required ",
           "(row names are not numeric).")
    }
  } else {
    ages_num <- as.numeric(ages)
    if (length(ages_num) != nrow(rate_mat_raw) || anyNA(ages_num))
      stop("`ages` must be a numeric vector of length nrow(rate_mat_raw).")
  }
  
  # --- years ---
  if (is.null(years)) {
    yrs_num <- suppressWarnings(as.numeric(colnames(rate_mat_raw)))
    if (anyNA(yrs_num) || length(yrs_num) == 0L) {
      stop("stmomo_forecast_iter(): numeric `years` vector required ",
           "(column names are not numeric).")
    }
  } else {
    yrs_num <- as.numeric(years)
    if (length(yrs_num) != ncol(rate_mat_raw) || anyNA(yrs_num))
      stop("`years` must be a numeric vector of length ncol(rate_mat_raw).")
  }
  yrs_w <- yrs_num[start_idx:end_idx]
  
  # --- Poisson counts / exposures ---
  Dxt <- round(rate_w * pop_w)
  Ext <- pop_w
  
  invalid <- !is.finite(Dxt) | !is.finite(Ext) | Ext <= 0 | Dxt < 0
  Dxt[invalid] <- NA_real_
  Ext[invalid] <- NA_real_
  
  wxt <- matrix(as.numeric(!invalid), nrow = nrow(Dxt), ncol = ncol(Dxt))
  dimnames(wxt) <- dimnames(Dxt)
  
  # Down-weight weakly-identified boundary cohorts for APC only.
  if (model_type == "APC") {
    wxt_cohort <- tryCatch(
      StMoMo::genWeightMat(ages = ages_num, years = yrs_w,
                           clip = as.integer(cohort_clip)),
      error = function(e) NULL
    )
    if (!is.null(wxt_cohort) &&
        all(dim(wxt_cohort) == dim(wxt))) {
      wxt <- wxt * wxt_cohort
    }
  }
  
  mod <- switch(model_type,
                "LC"  = StMoMo::lc(link  = "log"),
                "APC" = StMoMo::apc(link = "log"))
  
  # Fit using Dxt/Ext directly (the officially documented path).
  fit_res <- tryCatch(
    suppressWarnings(
      StMoMo::fit(mod,
                  Dxt     = Dxt,
                  Ext     = Ext,
                  ages    = ages_num,
                  years   = yrs_w,
                  wxt     = wxt,
                  verbose = verbose_fit)
    ),
    error = function(e) NULL
  )
  if (is.null(fit_res))
    return(rep(NA_real_, nrow(rate_mat_raw)))
  
  # Convergence / failure check. We only reject on an *explicit*
  # non-convergence or failure flag; a missing slot is treated as
  # "no information", not as failure, for robustness across StMoMo
  # versions.
  if (isTRUE(fit_res$fail) || identical(fit_res$conv, FALSE))
    return(rep(NA_real_, nrow(rate_mat_raw)))
  
  fc <- tryCatch(
    suppressWarnings(forecast(fit_res, h = h)),
    error = function(e) {
      message("  [", model_type, "] forecast() failed at h=", h,
              " (end_idx=", end_idx, "): ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fc))
    return(rep(NA_real_, nrow(rate_mat_raw)))
  
  # ---- robust extraction of h-step forecast rates ----
  # forecast.fitStMoMo returns $rates as an (ages x h) matrix in the
  # usual case, but `||` in R does NOT short-circuit on NA, and certain
  # StMoMo/forecast configurations can produce $rates with an NA dim,
  # a 3D (ages x h x sims) array, or a plain numeric vector when
  # h == 1. Handle each case explicitly.
  rates_obj <- fc$rates
  if (is.null(rates_obj)) rates_obj <- fc[["rate"]]  # defensive alias
  if (is.null(rates_obj)) {
    message("  [", model_type, "] forecast() returned no $rates at h=", h,
            "; names(fc) = ", paste(names(fc), collapse = ","))
    return(rep(NA_real_, nrow(rate_mat_raw)))
  }
  
  rates_mat <- if (is.matrix(rates_obj)) {
    rates_obj
  } else if (is.array(rates_obj) && length(dim(rates_obj)) == 3L) {
    apply(rates_obj, c(1L, 2L), mean, na.rm = TRUE)      # collapse sims
  } else if (is.numeric(rates_obj) && is.null(dim(rates_obj))) {
    matrix(rates_obj,
           nrow = nrow(rate_mat_raw),
           ncol = max(1L, length(rates_obj) %/% nrow(rate_mat_raw)))
  } else {
    NULL
  }
  
  nc <- if (is.matrix(rates_mat)) ncol(rates_mat) else NA_integer_
  if (is.null(rates_mat) || !is.finite(nc) || nc < h) {
    message("  [", model_type, "] unusable $rates at h=", h,
            " (class=", paste(class(rates_obj), collapse = "/"),
            ", dim=", paste(dim(rates_obj), collapse = "x"),
            ", nc=", nc, ")")
    return(rep(NA_real_, nrow(rate_mat_raw)))
  }
  
  rates_h <- rates_mat[, h]
  if (isTRUE(log_output)) {
    out <- suppressWarnings(log(rates_h))
  } else {
    out <- rates_h
  }
  out[!is.finite(out)] <- NA_real_
  out
}

# Convenience wrappers (names parallel the FTS code). Both resolve
# ages / years from dimnames (with sequential fallback + warning) so
# callers can pass a bare rate matrix without threading labels through.
stmomo_rolling_forecast_iteration <- function(rate_mat_raw, pop_mat,
                                              start_idx, end_idx, h,
                                              model_type,
                                              log_output  = TRUE,
                                              ages        = NULL,
                                              years       = NULL,
                                              cohort_clip = 3L) {
  ay <- resolve_ages_years_stmomo(rate_mat_raw, ages = ages, years = years)
  stmomo_forecast_iter(rate_mat_raw, pop_mat, start_idx, end_idx, h,
                       model_type  = model_type,
                       log_output  = log_output,
                       ages        = ay$ages,
                       years       = ay$years,
                       cohort_clip = cohort_clip)
}
stmomo_expanding_forecast_iteration <- function(rate_mat_raw, pop_mat,
                                                end_idx, h,
                                                model_type,
                                                log_output  = TRUE,
                                                ages        = NULL,
                                                years       = NULL,
                                                cohort_clip = 3L) {
  ay <- resolve_ages_years_stmomo(rate_mat_raw, ages = ages, years = years)
  stmomo_forecast_iter(rate_mat_raw, pop_mat, 1L, end_idx, h,
                       model_type  = model_type,
                       log_output  = log_output,
                       ages        = ay$ages,
                       years       = ay$years,
                       cohort_clip = cohort_clip)
}

# =====================================================================
# 2) POINT-FORECAST EVALUATION ON THE TEST SET
# =====================================================================
compare_test_performance_point_stmomo <- function(forecast_list,
                                                  actual_test_data,
                                                  apply_exp_transform) {
  out <- data.frame()
  scale_msg <- if (apply_exp_transform) "RAW (exp/additive)" else "native (additive)"
  cat(sprintf("Comparing point performance on test set (%s)...\n", scale_msg))
  for (method in names(forecast_list)) {
    pf  <- forecast_list[[method]]
    ok  <- is.finite(actual_test_data) & is.finite(pf$mean)
    y   <- as.vector(actual_test_data[ok])
    yh  <- as.vector(pf$mean[ok])
    rmse <- if (length(y) > 0) sqrt(mean((y - yh)^2)) else NA_real_
    mae  <- if (length(y) > 0) mean(abs(y - yh))      else NA_real_
    out <- rbind(out, data.frame(
      Method = method, RMSE = rmse, MAE = mae,
      Exp_Transform = apply_exp_transform,
      row.names = NULL
    ))
  }
  out
}

# =====================================================================
# 3) MAIN WRAPPER  -  one model, one horizon h
# =====================================================================
# Arguments:
#   rate_mat_raw : raw rate matrix, rates m_{x,t} or f_{x,t}  [age x year]
#                  (StMoMo fitting is always performed on the raw rate
#                   scale via Poisson counts D = round(rate * E).)
#   pop_mat      : central exposures E_{x,t}                   [age x year]
#   data_smooth  : actuals on the EVALUATION scale (the same scale the
#                  existing FTS code evaluates on):
#                    mortality -> log(smoothed m_{x,t})
#                    ASFR      -> smoothed ASFR (no log)
#   data_raw     : actuals on the EVALUATION scale for the raw series:
#                    mortality -> log(raw m_{x,t})
#                    ASFR      -> raw ASFR (no log)
#   log_output   : TRUE  => forecasts returned on log-rate scale (mortality)
#                  FALSE => forecasts returned on raw-rate scale  (ASFR)
#   mode         : "expanding" (default) or "rolling"
#   apply_exp_transform :
#                  TRUE  => evaluate on exp(eval-scale); FALSE => on eval-scale.
#                  Requires log_output = TRUE (enforced).
#                  Matches the FTS function semantics. In line with the
#                  existing mortality / ASFR pipelines, FALSE is recommended.
#   ages, years  : OPTIONAL numeric vectors. Strongly recommended for ASFR
#                  because the APC cohort index c = t - x depends on the
#                  actual age origin (e.g., 15-49), not on the row index.
#   cohort_clip  : integer passed to StMoMo::genWeightMat for APC only.
# =====================================================================
forecast_stmomo_point <- function(rate_mat_raw, pop_mat,
                                  data_smooth,
                                  data_raw,
                                  end_train_index,
                                  h,
                                  model_type = c("LC", "APC"),
                                  mode       = c("expanding", "rolling"),
                                  log_output = TRUE,
                                  apply_exp_transform = FALSE,
                                  ages        = NULL,
                                  years       = NULL,
                                  cohort_clip = 3L,
                                  seed = 42) {
  model_type <- match.arg(model_type)
  mode       <- match.arg(mode)
  
  # Sanity guard: exp() only makes sense if forecasts are on log scale.
  if (isTRUE(apply_exp_transform) && !isTRUE(log_output)) {
    stop("apply_exp_transform = TRUE requires log_output = TRUE ",
         "(cannot exponentiate rate-scale forecasts).")
  }
  
  # Resolve ages/years once at the wrapper level. This is the same logic
  # the driver scripts use; bringing it inside the wrapper removes the
  # requirement that every caller thread them through manually.
  ay       <- resolve_ages_years_stmomo(rate_mat_raw, ages = ages, years = years)
  ages_r   <- ay$ages
  years_r  <- ay$years
  
  set.seed(seed)
  
  num_rows    <- nrow(data_smooth)
  num_periods <- ncol(data_smooth)
  train_len   <- end_train_index
  
  # ------------- split (identical to Point_forecast_function.R) -------------
  available_all <- (end_train_index + h):num_periods
  if (length(available_all) < 2L) stop("Not enough data for validation and test sets.")
  
  total_start <- num_periods - end_train_index
  base_pairs  <- floor(total_start / 2L)
  pairs_h     <- max(base_pairs - (h - 1L), 1L)
  val_len  <- pairs_h
  test_len <- pairs_h
  
  val_pos        <- seq_len(val_len)
  test_start_pos <- val_len + h
  test_pos       <- seq.int(test_start_pos, test_start_pos + test_len - 1L)
  if (max(test_pos) > length(available_all))
    stop(sprintf("Split overflow: need %d positions, but only %d available (h=%d).",
                 max(test_pos), length(available_all), h))
  
  validation_target_indices <- available_all[val_pos]
  test_target_indices       <- available_all[test_pos]
  
  cat(sprintf(
    "[%s / %s] h=%d | total_start=%d | val=%d (%d..%d) | test=%d (%d..%d) | gap=%d\n",
    model_type, mode, h, total_start,
    val_len,  min(validation_target_indices), max(validation_target_indices),
    test_len, min(test_target_indices),       max(test_target_indices),
    h - 1L))
  
  # ------------- test forecasts -------------
  n_test <- length(test_target_indices)
  fc_test <- matrix(NA_real_, nrow = num_rows, ncol = n_test)
  for (i in 1:n_test) {
    end_idx   <- test_target_indices[i] - h
    start_idx <- if (mode == "rolling") max(1L, end_idx - train_len + 1L) else 1L
    fc_test[, i] <- stmomo_forecast_iter(
      rate_mat_raw, pop_mat, start_idx, end_idx, h,
      model_type  = model_type,
      log_output  = log_output,
      ages        = ages_r,
      years       = years_r,
      cohort_clip = cohort_clip
    )
  }
  
  test_data_raw_actuals    <- data_raw[,    test_target_indices, drop = FALSE]
  test_data_smooth_actuals <- data_smooth[, test_target_indices, drop = FALSE]
  
  if (isTRUE(apply_exp_transform)) {
    fc_test_eval     <- exp(fc_test)
    actuals_for_eval <- exp(test_data_raw_actuals)
  } else {
    fc_test_eval     <- fc_test
    actuals_for_eval <- test_data_smooth_actuals
  }
  
  test_forecast_list <- list()
  test_forecast_list[[model_type]] <- list(mean = fc_test_eval)
  
  cat(sprintf("\nEvaluating %s on %s scale...\n", model_type,
              if (apply_exp_transform) "RAW (exp/additive)" else "native (additive)"))
  test_results <- compare_test_performance_point_stmomo(
    forecast_list       = test_forecast_list,
    actual_test_data    = actuals_for_eval,
    apply_exp_transform = apply_exp_transform
  )
  test_results$Model <- model_type
  test_results$Mode  <- mode
  
  cat("--- forecast_stmomo_point finished (", model_type, ", ", mode, ") ---\n", sep = "")
  list(
    settings = list(h = h, model_type = model_type, mode = mode,
                    log_output = log_output,
                    exp_transform_applied = apply_exp_transform),
    test_phase_performance = test_results,
    test_phase_forecasts = list(
      actuals_raw    = test_data_raw_actuals,
      actuals_smooth = test_data_smooth_actuals,
      forecast       = fc_test       # same scale as log_output flag
    )
  )
}

# =====================================================================
# 4) EXAMPLE USAGE (POINT FORECAST)  ---- not run ----
# =====================================================================
if (TRUE) {
  # ---- Mortality example ----
  demog_list <- readRDS("demog_list.rds")
  smooth_mort_f <- demog_list$AUS$Female_smooth$rate[[1]]
  raw_mort_f    <- demog_list$AUS$Female_raw$rate[[1]]
  pop_f         <- demog_list$AUS$Female_raw$pop[[1]]
  
  res_lc_mort <- forecast_stmomo_point(
    rate_mat_raw    = raw_mort_f,
    pop_mat         = pop_f,
    data_smooth     = log(smooth_mort_f),  # log scale (mortality)
    data_raw        = log(raw_mort_f),     # log scale (mortality)
    end_train_index = ncol(smooth_mort_f) - 40L,
    h               = 1,
    model_type      = "LC",
    mode            = "expanding",
    log_output      = TRUE,                # mortality -> log forecasts
    apply_exp_transform = FALSE
  )
  print(res_lc_mort$test_phase_performance)
  
  # ---- ASFR example ----
  demog_list_asfr <- readRDS("demog_list_asfr.rds")
  raw_asfr    <- demog_list_asfr$CAN$raw$rate[[1]]
  smooth_asfr <- demog_list_asfr$CAN$smooth$rate[[1]]
  exp_asfr    <- demog_list_asfr$CAN$raw$pop[[1]]
  
  res_apc_asfr <- forecast_stmomo_point(
    rate_mat_raw    = raw_asfr,
    pop_mat         = exp_asfr,
    data_smooth     = smooth_asfr,         # rate scale (no log for ASFR)
    data_raw        = raw_asfr,
    end_train_index = ncol(smooth_asfr) - 40L,
    h               = 1,
    model_type      = "APC",
    mode            = "expanding",
    log_output      = FALSE,               # ASFR -> rate-scale forecasts
    apply_exp_transform = FALSE
  )
  print(res_apc_asfr$test_phase_performance)
}