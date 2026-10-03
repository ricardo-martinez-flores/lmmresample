#' Calibrate a permutation test by simulation
#'
#' Checks whether [perm_test()] controls the type I error rate for the
#' user's own design and model. Data are simulated repeatedly under the null
#' hypothesis from the fitted model, with the tested coefficient set to zero,
#' and the permutation test is applied to each simulated data set. If the test
#' is well calibrated, the proportion of rejections is close to `alpha` and
#' the p-values are uniformly distributed. For comparison, the rejection rate
#' of the model's Wald test is reported for the same simulated data sets.
#'
#' @details
#' Simulated responses keep everything else from the fitted model: the
#' remaining fixed effects, the estimated random-effects covariance and the
#' residual variance. Two null hypotheses are available:
#'
#' * `null = "mean"` (default): the average effect of `term` is zero, but
#'   participants differ in their individual effect as estimated by the random
#'   slope of `term`. This is the hypothesis usually claimed in applied work,
#'   and the one under which a permutation test can fail, so it is the
#'   informative check.
#' * `null = "sharp"`: the effect of `term` is zero for every participant; the
#'   random slope of `term` is also set to zero. The permutation test is exact
#'   under this null when units are exchangeable, so it serves as a reference.
#'
#' Residuals are simulated as independent by default. For densely sampled
#' series, `ar1` adds first-order autocorrelation within each `series` (e.g.
#' trial), ordered by `time`, keeping the residual standard deviation of the
#' model. This reproduces the dependence that makes standard inference
#' unreliable for such data.
#'
#' The check is only as realistic as the fitted model: anything the model
#' leaves in the residuals is simulated as random noise. In particular, a
#' common time course that the model does not capture (e.g. the average
#' response curve, fitted only as a linear trend) is treated as noise that
#' varies from trial to trial, which can make the simulated data behave
#' quite differently from the real data. Fit the time course flexibly
#' (e.g. with splines or time as a factor) before calibrating; residual
#' time courses can be inspected with [diag_timecourse()].
#'
#' Each simulation runs a complete permutation test, so the total number of
#' model fits is `n_sim * (B + 1)`. Simulations run in parallel through the
#' \pkg{future} framework (see [perm_test()]). The defaults give a quick
#' check; for a definitive calibration increase `n_sim` (e.g. 1000), which
#' narrows the confidence interval of the rejection rate.
#'
#' @inheritParams perm_test
#' @param method Resampling method, as in [perm_test()].
#' @param null Null hypothesis under which data are simulated: `"mean"` or
#'   `"sharp"`. See Details.
#' @param ar1 Optional lag-1 autocorrelation of the simulated residuals within
#'   each `series`, in \eqn{[0, 1)}.
#' @param series,time Names of the columns identifying each series (e.g.
#'   `"trial"`) and the time order within it. Required when `ar1` is
#'   supplied; `time` is also used to match samples when Freedman-Lane
#'   residuals are permuted between units (see [perm_test()]).
#' @param n_sim Number of simulated data sets.
#' @param B Number of permutations per simulated data set.
#' @param alpha Significance level at which rejection rates are evaluated.
#'
#' @return An object of class `lmmr_calib` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()]. Its elements include `p_perm`
#'   and `p_wald` (p-values of the permutation and Wald tests in each
#'   simulation) and `rejection` (rejection rates at `alpha` with exact
#'   binomial confidence intervals).
#'
#' @seealso [perm_test()]
#'
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 10,
#'                 n_trials = 8, seed = 1)
#' m <- lm(y ~ condition + participant, data = d)
#'
#' # Small values to keep the example fast
#' cal <- perm_calibrate(m, term = "condition", unit = "trial",
#'                       exchange = exch_within("participant"),
#'                       n_sim = 20, B = 49, seed = 1)
#' cal
#' @export
perm_calibrate <- function(model,
                           term,
                           unit,
                           exchange,
                           coef = NULL,
                           method = c("auto", "relabel", "freedman-lane"),
                           null = c("mean", "sharp"),
                           ar1 = NULL,
                           series = NULL,
                           time = NULL,
                           n_sim = 200,
                           B = 199,
                           alpha = 0.05,
                           seed = NULL,
                           data = NULL) {
  call <- match.call()
  check_model(model)
  check_column_name(term, "term")
  method <- match.arg(method)
  null <- match.arg(null)
  check_count(n_sim, "n_sim", min = 1)
  check_count(B, "B", min = 1)
  check_number(alpha, "alpha", min = 0, max = 1, min_inclusive = FALSE)
  if (!is.null(ar1)) {
    if (is.list(ar1) && !is.null(ar1$phi)) ar1 <- ar1$phi
    check_number(ar1, "ar1", min = 0)
    if (ar1 >= 1) cli::cli_abort("{.arg ar1} must be smaller than 1.")
  }

  data <- model_data(model, data)
  test <- resolve_test(model, data, term, coef, method)
  engine <- build_engine(model, data, term, test, unit, exchange, time)
  check_permutation_space(engine$space, B, exchange, unit)

  response <- response_name(model)
  simulate_y <- make_null_simulator(model, data, test$coefs, null, ar1,
                                    series, time)
  refit <- make_refitter(model)
  label <- stat_label(model)
  df_resid <- if (inherits(model, "lm")) stats::df.residual(model) else Inf

  sim_seeds <- with_seed(seed, sample.int(.Machine$integer.max, n_sim))
  rng_kind <- RNGkind()

  use_progress <- requireNamespace("progressr", quietly = TRUE)
  if (use_progress) p <- progressr::progressor(steps = n_sim)
  results <- future.apply::future_lapply(seq_len(n_sim), function(i) {
    out <- with_seed(sim_seeds[i], kind = rng_kind, expr = {
      sim_data <- data
      sim_data[[response]] <- simulate_y()
      calibrate_one(model, refit, sim_data, term, test, unit, exchange, time, B,
                    df_resid, label)
    })
    if (use_progress) p()
    out
  }, future.seed = FALSE, future.packages = refit_packages(model))

  p_perm <- vapply(results, function(r) r$p_perm, numeric(1))
  p_wald <- vapply(results, function(r) r$p_wald, numeric(1))
  b_used <- vapply(results, function(r) r$B_used, numeric(1))

  structure(
    list(
      p_perm = p_perm,
      p_wald = p_wald,
      B_used = b_used,
      rejection = rejection_table(list(Permutation = p_perm, Wald = p_wald),
                                  alpha),
      uniformity = uniformity_test(p_perm),
      alpha = alpha,
      null = null,
      ar1 = ar1,
      series = series,
      term = term,
      coef = test$coefs,
      method = test$method,
      unit = unit,
      exchange = exchange,
      n_sim = n_sim,
      B = B,
      seed = seed,
      call = call
    ),
    class = "lmmr_calib"
  )
}

# One simulated data set: observed statistic, Wald p-value and permutation
# p-value, with refits run sequentially (parallelism is over simulations).
calibrate_one <- function(model, refit, sim_data, term, test, unit, exchange,
                          time, B, df_resid, label) {
  failed <- list(p_perm = NA_real_, p_wald = NA_real_, B_used = 0)
  obs <- safe_refit(refit, sim_data, test$coefs, test$type)
  if (obs$status %in% c("failed", "nonconverged")) return(failed)
  stat_obs <- obs$stat[[1]]
  p_wald <- if (test$type == "chi2") {
    stats::pchisq(stat_obs, df = length(test$coefs), lower.tail = FALSE)
  } else if (label == "t" && is.finite(df_resid)) {
    2 * stats::pt(-abs(stat_obs), df_resid)
  } else {
    2 * stats::pnorm(-abs(stat_obs))
  }
  engine <- tryCatch(
    suppressWarnings(build_engine(model, sim_data, term, test, unit,
                                  exchange, time)),
    error = function(e) NULL
  )
  if (is.null(engine)) return(failed)
  draws <- engine$draw(B)
  null <- vapply(seq_len(B), function(b) {
    res <- safe_refit(refit, engine$make_data(draws, b), test$coefs,
                      test$type)
    if (res$status %in% c("failed", "nonconverged")) NA_real_ else res$stat[[1]]
  }, numeric(1))
  null <- null[!is.na(null)]
  alt <- if (test$type == "chi2") "greater" else "two.sided"
  list(p_perm = perm_pvalue(stat_obs, null, alt), p_wald = p_wald,
       B_used = length(null))
}

rejection_table <- function(pvalues, alpha) {
  rows <- lapply(names(pvalues), function(method) {
    p <- pvalues[[method]]
    p <- p[!is.na(p)]
    k <- sum(p <= alpha)
    ci <- stats::binom.test(k, length(p))$conf.int
    data.frame(method = method, rejection = k / length(p), lower = ci[1],
               upper = ci[2], n = length(p))
  })
  do.call(rbind, rows)
}

uniformity_test <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) < 5) return(NA_real_)
  suppressWarnings(stats::ks.test(p, "punif")$p.value)
}

calib_verdict <- function(row, alpha) {
  if (row$lower <= alpha && row$upper >= alpha) return("consistent with alpha")
  if (row$lower > alpha) "liberal (too many rejections)" else
    "conservative (too few rejections)"
}

# Methods ---------------------------------------------------------------------

#' @export
print.lmmr_calib <- function(x, ...) {
  cat("\nCalibration of the permutation test for `", x$term, "`\n\n", sep = "")
  null_text <- if (x$null == "mean") {
    "mean effect zero, participant variation kept"
  } else {
    "no effect in any participant"
  }
  cat("Null:        ", null_text, "\n", sep = "")
  cat("Residuals:   ", if (is.null(x$ar1)) "independent" else
        paste0("AR(1), phi = ", format(round(x$ar1, 3), nsmall = 3),
               " within `", x$series, "`"),
      "\n", sep = "")
  cat("Exchange:    ", format(x$exchange, unit = x$unit), "\n", sep = "")
  cat("Simulations: ", x$n_sim, ", with ", x$B, " permutations each\n\n",
      sep = "")
  cat("Rejection rate at alpha = ", x$alpha, ":\n", sep = "")
  for (i in seq_len(nrow(x$rejection))) {
    row <- x$rejection[i, ]
    cat(sprintf("  %-12s %.3f  [%.3f, %.3f]  %s\n", row$method, row$rejection,
                row$lower, row$upper, calib_verdict(row, x$alpha)))
  }
  if (!is.na(x$uniformity)) {
    cat("\nUniformity of permutation p-values (Kolmogorov-Smirnov): p ",
        format_p(x$uniformity, prefix = TRUE), "\n", sep = "")
  }
  n_failed <- sum(is.na(x$p_perm))
  if (n_failed > 0) {
    cat(n_failed, " simulation", if (n_failed > 1) "s", " failed to fit\n",
        sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Tidy a calibration
#'
#' @param x An object of class `lmmr_calib`.
#' @param type `"summary"` returns the rejection rates; `"pvalues"` returns
#'   the p-values of each simulation.
#' @param ... Not used.
#' @return A data frame.
#' @export
tidy.lmmr_calib <- function(x, type = c("summary", "pvalues"), ...) {
  type <- match.arg(type)
  if (type == "pvalues") {
    return(data.frame(simulation = seq_along(x$p_perm), p_perm = x$p_perm,
                      p_wald = x$p_wald))
  }
  x$rejection
}

#' Plot a calibration
#'
#' @param x An object of class `lmmr_calib`.
#' @param type `"pp"` compares the observed p-values with the uniform
#'   distribution expected under a well-calibrated test, with a pointwise 95%
#'   band. `"rejection"` shows the rejection rate of each method with its
#'   confidence interval.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_calib <- function(x, type = c("pp", "rejection"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  type <- match.arg(type)
  colours <- c(Permutation = "grey15", Wald = "#C0392B")

  if (type == "rejection") {
    rej <- x$rejection
    rej$method <- factor(rej$method, levels = names(colours))
    return(
      ggplot2::ggplot(rej, ggplot2::aes(x = .data$method, y = .data$rejection,
                                        colour = .data$method)) +
        ggplot2::geom_hline(yintercept = x$alpha, linetype = "dashed",
                            colour = "grey40") +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lower,
                                              ymax = .data$upper),
                                 linewidth = 0.8, size = 0.6) +
        ggplot2::scale_colour_manual(values = colours, guide = "none") +
        ggplot2::scale_y_continuous(limits = c(0, NA)) +
        ggplot2::labs(
          x = NULL, y = "Rejection rate under the null",
          title = paste0("Type I error: ", x$term),
          subtitle = paste0("Dashed line: alpha = ", x$alpha, "; ", x$n_sim,
                            " simulations, 95% exact intervals")
        ) +
        theme_lmmr()
    )
  }

  pp <- do.call(rbind, lapply(c("Permutation", "Wald"), function(method) {
    p <- sort(if (method == "Permutation") x$p_perm else x$p_wald)
    n <- length(p)
    data.frame(method = method, expected = seq_len(n) / (n + 1), observed = p)
  }))
  pp$method <- factor(pp$method, levels = names(colours))
  n <- sum(!is.na(x$p_perm))
  i <- seq_len(n)
  band <- data.frame(expected = i / (n + 1),
                     lower = stats::qbeta(0.025, i, n - i + 1),
                     upper = stats::qbeta(0.975, i, n - i + 1))
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band,
                         ggplot2::aes(x = .data$expected, ymin = .data$lower,
                                      ymax = .data$upper),
                         fill = "grey88") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50",
                         linetype = "dashed") +
    ggplot2::geom_step(data = pp,
                       ggplot2::aes(x = .data$expected, y = .data$observed,
                                    colour = .data$method),
                       linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Expected p-value (uniform)", y = "Observed p-value",
      title = paste0("Calibration of p-values: ", x$term),
      subtitle = "A calibrated test follows the diagonal (band: 95% pointwise)"
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
