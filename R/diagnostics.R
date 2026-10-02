#' Autocorrelation of residuals within series
#'
#' Estimates the autocorrelation of the model residuals within each series
#' (e.g. the samples of a trial), ordered by time, and averages it across
#' series. Strong autocorrelation is the reason why Wald inference on densely
#' sampled data is unreliable and why whole series should be resampled. The
#' lag-1 estimate can be passed to [perm_calibrate()] to simulate data with
#' the same dependence.
#'
#' @details
#' Each series is centred on its own mean before computing its
#' autocorrelation, which biases estimates slightly towards zero for short
#' series (by roughly \eqn{(1 + 4\phi) / n} at lag 1 for a series of length
#' \eqn{n}). Series shorter than `lag_max + 2` samples are skipped.
#'
#' @param model A linear model fitted with [stats::lm()] or [lme4::lmer()].
#' @param series Name of the column identifying each series (e.g. `"trial"`).
#' @param time Name of the column giving the order of samples within series.
#' @param lag_max Maximum lag. Defaults to 20 or the length of the shortest
#'   series minus 2, whichever is smaller.
#' @param data The data used to fit `model`, if it cannot be recovered from
#'   the model call.
#'
#' @return An object of class `lmmr_acf` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()]. Its element `phi` is the mean
#'   lag-1 autocorrelation, and `acf` a data frame with the mean, 10th and
#'   90th percentiles of the autocorrelation across series at each lag.
#'
#' @seealso [perm_calibrate()], [diag_timecourse()]
#'
#' @examples
#' d <- sim_blocks("within", n_participants = 6, n_trials = 10, n_time = 40,
#'                 ar1 = 0.8, seed = 1)
#' m <- lm(y ~ condition + participant, data = d)
#' ac <- diag_acf(m, series = "trial", time = "time")
#' ac
#' @export
diag_acf <- function(model, series, time, lag_max = NULL, data = NULL) {
  check_model(model)
  check_column_name(series, "series")
  check_column_name(time, "time")
  data <- model_data(model, data)
  missing <- setdiff(c(series, time), names(data))
  if (length(missing) > 0) {
    cli::cli_abort("Column{?s} {.field {missing}} not found in the data.")
  }
  check_series_ids(data, series, time)
  r <- as.vector(stats::residuals(model))
  ord <- order(data[[series]], data[[time]])
  groups <- split(r[ord], as.character(data[[series]])[ord])
  lengths <- lengths(groups)
  if (is.null(lag_max)) lag_max <- max(1, min(20, min(lengths) - 2))
  check_count(lag_max, "lag_max", min = 1)
  usable <- groups[lengths >= lag_max + 2]
  if (length(usable) == 0) {
    cli::cli_abort("No series has at least {lag_max + 2} samples; reduce
                    {.arg lag_max}.")
  }
  acfs <- vapply(usable, function(x) {
    if (stats::sd(x) == 0) return(rep(NA_real_, lag_max))
    as.vector(stats::acf(x, lag.max = lag_max, plot = FALSE,
                         demean = TRUE)$acf)[-1]
  }, numeric(lag_max))
  acfs <- matrix(acfs, nrow = lag_max)
  tab <- data.frame(
    lag = seq_len(lag_max),
    acf = rowMeans(acfs, na.rm = TRUE),
    lower = apply(acfs, 1, stats::quantile, 0.1, na.rm = TRUE),
    upper = apply(acfs, 1, stats::quantile, 0.9, na.rm = TRUE)
  )
  structure(
    list(acf = tab, phi = tab$acf[1], series = series, time = time,
         n_series = length(usable), n_skipped = length(groups) - length(usable),
         median_length = stats::median(lengths)),
    class = "lmmr_acf"
  )
}

#' @export
print.lmmr_acf <- function(x, ...) {
  cat("\nResidual autocorrelation within `", x$series, "` (ordered by `",
      x$time, "`)\n\n", sep = "")
  cat("Series:            ", x$n_series, " (median length ",
      x$median_length, ")", if (x$n_skipped > 0) {
        paste0("; ", x$n_skipped, " too short and skipped")
      }, "\n", sep = "")
  cat("Lag-1 correlation: ", formatC(x$phi, digits = 3, format = "f"),
      "\n", sep = "")
  show <- x$acf[x$acf$lag %in% c(1, 2, 5, 10, 20), ]
  cat("Mean autocorrelation at lags ",
      paste(show$lag, collapse = ", "), ": ",
      paste(formatC(show$acf, digits = 2, format = "f"), collapse = ", "),
      "\n", sep = "")
  if (x$phi > 0.3) {
    cat("\nResiduals are strongly autocorrelated: Wald tests that treat",
        "samples as\nindependent are unreliable. Resample whole series",
        "(see perm_test()).\n")
  }
  cat("\n")
  invisible(x)
}

#' Tidy a residual autocorrelation diagnostic
#' @param x An object of class `lmmr_acf`.
#' @param ... Not used.
#' @return A data frame with one row per lag.
#' @export
tidy.lmmr_acf <- function(x, ...) x$acf

#' Plot a residual autocorrelation diagnostic
#' @param x An object of class `lmmr_acf`.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_acf <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  ggplot2::ggplot(x$acf, ggplot2::aes(x = .data$lag, y = .data$acf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_line(colour = "grey15", linewidth = 0.7) +
    ggplot2::geom_point(colour = "grey15", size = 1.5) +
    ggplot2::scale_y_continuous(limits = c(min(-0.1, x$acf$lower), 1)) +
    ggplot2::labs(
      x = "Lag (samples)", y = "Autocorrelation of residuals",
      title = paste0("Residual autocorrelation within ", x$series),
      subtitle = paste0("Mean across ", x$n_series,
                        " series; band: 10th-90th percentile. Lag-1 = ",
                        round(x$phi, 2))
    ) +
    theme_lmmr()
}

#' Residual time courses by condition
#'
#' Averages the model residuals over time, separately for each level of a
#' grouping variable such as the condition. If the model describes the time
#' course well, the mean residual stays close to zero at all times in every
#' condition. Systematic departures indicate that the shape of the response
#' (e.g. a linear term for a curved response) or its difference between
#' conditions (a missing condition by time interaction) is misspecified.
#'
#' @details
#' A misspecified time course does not invalidate permutation tests of the
#' model terms, but the coefficients may then not answer the intended
#' question and the tests lose power. Flexible terms such as
#' `splines::ns(time, df)` and their interaction with condition are usually
#' adequate. The mean residual and its standard error are computed across
#' clusters (e.g. participants) when `cluster` is given, and across rows
#' otherwise.
#'
#' @inheritParams diag_acf
#' @param by Optional name of a grouping column (e.g. `"condition"`).
#' @param cluster Optional name of the column identifying independent
#'   clusters (e.g. `"participant"`), used for the standard errors.
#' @param bins Number of time bins used when `time` has more than 50 distinct
#'   values. `NULL` keeps the original time points.
#'
#' @return An object of class `lmmr_timecourse` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()].
#' @seealso [diag_acf()]
#'
#' @examples
#' d <- sim_blocks("within", n_participants = 8, n_trials = 10, n_time = 30,
#'                 effect_condition = 1, effect_onset = 0.3, seed = 1)
#' m <- lm(y ~ condition + time + participant, data = d)
#' tc <- diag_timecourse(m, time = "time", by = "condition",
#'                       cluster = "participant")
#' tc
#' @export
diag_timecourse <- function(model, time, by = NULL, cluster = NULL,
                            bins = 30, data = NULL) {
  check_model(model)
  check_column_name(time, "time")
  if (!is.null(by)) check_column_name(by, "by")
  if (!is.null(cluster)) check_column_name(cluster, "cluster")
  data <- model_data(model, data)
  missing <- setdiff(c(time, by, cluster), names(data))
  if (length(missing) > 0) {
    cli::cli_abort("Column{?s} {.field {missing}} not found in the data.")
  }
  r <- as.vector(stats::residuals(model))
  t <- data[[time]]
  if (!is.numeric(t)) cli::cli_abort("{.field {time}} must be numeric.")
  if (!is.null(bins) && length(unique(t)) > 50) {
    check_count(bins, "bins", min = 2)
    breaks <- seq(min(t), max(t), length.out = bins + 1)
    mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
    t <- mids[findInterval(t, breaks, rightmost.closed = TRUE,
                           all.inside = TRUE)]
  }
  grp <- if (is.null(by)) rep("all", length(r)) else as.character(data[[by]])
  cl <- if (is.null(cluster)) seq_along(r) else as.character(data[[cluster]])

  # Mean per cluster first, then across clusters
  per_cluster <- stats::aggregate(r ~ cl + grp + t,
                                  data = data.frame(r = r, cl = cl, grp = grp,
                                                    t = t),
                                  FUN = mean)
  summ <- do.call(rbind, lapply(
    split(per_cluster, list(per_cluster$grp, per_cluster$t), drop = TRUE),
    function(x) {
      n <- nrow(x)
      data.frame(group = x$grp[1], time = x$t[1], mean = mean(x$r),
                 se = if (n > 1) stats::sd(x$r) / sqrt(n) else NA_real_, n = n)
    }
  ))
  summ <- summ[order(summ$group, summ$time), ]
  rownames(summ) <- NULL
  sd_r <- stats::sd(r)
  structure(
    list(timecourse = summ, time = time, by = by, cluster = cluster,
         residual_sd = sd_r,
         max_abs = tapply(abs(summ$mean), summ$group, max) / sd_r),
    class = "lmmr_timecourse"
  )
}

#' @export
print.lmmr_timecourse <- function(x, ...) {
  cat("\nResidual time course over `", x$time, "`",
      if (!is.null(x$by)) paste0(" by `", x$by, "`"), "\n\n", sep = "")
  cat("Largest mean residual, in residual SD units:\n")
  for (g in names(x$max_abs)) {
    cat(sprintf("  %-12s %.2f\n", g, x$max_abs[[g]]))
  }
  if (any(x$max_abs > 0.1)) {
    cat("\nThe mean residual departs from zero over time: consider flexible",
        "time terms\n(e.g. splines::ns(time, 4)) and their interaction with",
        "the condition.\n")
  }
  cat("\n")
  invisible(x)
}

#' Tidy a residual time-course diagnostic
#' @param x An object of class `lmmr_timecourse`.
#' @param ... Not used.
#' @return A data frame with the mean residual and its standard error at each
#'   time point and group.
#' @export
tidy.lmmr_timecourse <- function(x, ...) x$timecourse

#' Plot a residual time-course diagnostic
#' @param x An object of class `lmmr_timecourse`.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_timecourse <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  tc <- x$timecourse
  tc$lower <- tc$mean - 1.96 * tc$se
  tc$upper <- tc$mean + 1.96 * tc$se
  multi <- length(unique(tc$group)) > 1
  p <- ggplot2::ggplot(tc, ggplot2::aes(x = .data$time, y = .data$mean)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed")
  if (multi) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                        fill = .data$group), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(colour = .data$group), linewidth = 0.8) +
      ggplot2::scale_colour_manual(values = c("grey15", "#C0392B", "#2471A3",
                                              "#7D3C98", "#B7950B"),
                                   name = x$by) +
      ggplot2::scale_fill_manual(values = c("grey15", "#C0392B", "#2471A3",
                                            "#7D3C98", "#B7950B"),
                                 name = x$by)
  } else {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                           fill = "grey85") +
      ggplot2::geom_line(colour = "grey15", linewidth = 0.8)
  }
  p +
    ggplot2::labs(
      x = x$time, y = "Mean residual",
      title = "Residual time course",
      subtitle = paste0("Should stay near zero; band: 95% interval",
                        if (!is.null(x$cluster)) {
                          paste0(" across ", x$cluster, "s")
                        })
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
