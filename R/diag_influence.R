#' Leave-one-cluster-out influence
#'
#' Refits the model leaving out each cluster (e.g. participant) in turn and
#' reports how much each fixed-effect estimate changes. Large changes flag
#' influential participants whose removal would alter the conclusions.
#'
#' @details
#' The change is also expressed in units of the standard error of the full
#' model estimate. As a rough guide, values above 1 deserve inspection.
#' Refits run in parallel through the \pkg{future} framework.
#'
#' @param model A linear model fitted with [stats::lm()] or [lme4::lmer()].
#' @param cluster Name of the column identifying the clusters to leave out.
#' @param terms Names of the fixed-effect coefficients to report. By default
#'   all except the intercept.
#' @param data The data used to fit `model`, if it cannot be recovered from
#'   the model call.
#'
#' @return An object of class `lmmr_loso` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()]. Its element `estimates` is a
#'   matrix with one row per left-out cluster and one column per coefficient.
#'
#' @seealso [boot_ci()]
#'
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 10,
#'                 n_trials = 10, effect_condition = 0.5, seed = 1)
#' m <- lm(y ~ condition + participant, data = d)
#' loso <- diag_loso(m, cluster = "participant", terms = "conditionB")
#' loso
#' @export
diag_loso <- function(model, cluster, terms = NULL, data = NULL) {
  check_model(model)
  check_column_name(cluster, "cluster")
  data <- model_data(model, data)
  if (!cluster %in% names(data)) {
    cli::cli_abort("Column {.field {cluster}} not found in the data.")
  }
  est_full <- fixed_estimates(model)
  se_full <- sqrt(diag(as.matrix(stats::vcov(model))))[names(est_full)]
  terms <- check_terms(terms, est_full)

  refit <- make_refitter(model)
  res <- loso_refits(refit, data, cluster, terms)
  change <- sweep(res$estimates, 2, est_full[terms])
  structure(
    list(
      estimates = res$estimates,
      status = res$status,
      full = est_full[terms],
      se = se_full[terms],
      change = change,
      std_change = sweep(change, 2, se_full[terms], "/"),
      cluster = cluster
    ),
    class = "lmmr_loso"
  )
}

# Leave-one-cluster-out refits (also used for BCa acceleration).
loso_refits <- function(refit, data, cluster, coefs) {
  ids <- unique(as.character(data[[cluster]]))
  if (length(ids) < 3) {
    cli::cli_abort("At least 3 levels of {.field {cluster}} are needed.",
                   call = rlang::caller_env())
  }
  out <- future.apply::future_lapply(ids, function(id) {
    newdata <- droplevels(data[as.character(data[[cluster]]) != id, ,
                               drop = FALSE])
    safe_refit_estimates(refit, newdata, coefs)
  }, future.seed = FALSE)
  est <- matrix(vapply(out, function(o) o$est, numeric(length(coefs))),
                ncol = length(coefs), byrow = TRUE,
                dimnames = list(ids, coefs))
  status <- stats::setNames(vapply(out, function(o) o$status, character(1)),
                            ids)
  list(estimates = est, status = status)
}

# Refit and return fixed-effect estimates, recording problems.
safe_refit_estimates <- function(refit, newdata, coefs) {
  res <- safe_fit(refit, newdata)
  if (is.null(res$fit)) {
    return(list(est = stats::setNames(rep(NA_real_, length(coefs)), coefs),
                status = "failed"))
  }
  est <- fixed_estimates(res$fit)[coefs]
  names(est) <- coefs
  if (anyNA(est) && res$status == "ok") res$status <- "failed"
  if (res$status == "nonconverged") est[] <- NA_real_
  list(est = est, status = res$status)
}

fixed_estimates <- function(model) {
  if (is_mixed(model)) lme4::fixef(model) else stats::coef(model)
}

check_terms <- function(terms, est, call = rlang::caller_env()) {
  if (is.null(terms)) {
    terms <- setdiff(names(est)[!is.na(est)], "(Intercept)")
  }
  missing <- setdiff(terms, names(est))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "{.val {missing}} {?is/are} not {?a /}fixed-effect coefficient{?s}.",
      "i" = "Available: {.val {names(est)}}."
    ), call = call)
  }
  terms
}

#' @export
print.lmmr_loso <- function(x, digits = 3, ...) {
  cat("\nLeave-one-", x$cluster, "-out influence (", nrow(x$estimates),
      " refits)\n\n", sep = "")
  rows <- lapply(colnames(x$estimates), function(cf) {
    sc <- x$std_change[, cf]
    i <- which.max(abs(sc))
    data.frame(
      Coefficient = cf,
      Estimate = formatC(x$full[[cf]], digits = digits, format = "f"),
      Range = paste0(formatC(min(x$estimates[, cf], na.rm = TRUE),
                             digits = digits, format = "f"), " to ",
                     formatC(max(x$estimates[, cf], na.rm = TRUE),
                             digits = digits, format = "f")),
      `Most influential` = paste0(names(sc)[i], " (", sprintf("%.2f", sc[i]),
                                  " SE)"),
      check.names = FALSE
    )
  })
  print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
  n_bad <- sum(x$status %in% c("failed", "nonconverged"))
  if (n_bad > 0) cat("\n", n_bad, " refit(s) failed or did not converge\n",
                     sep = "")
  if (any(abs(x$std_change) > 1, na.rm = TRUE)) {
    cat("\nSome ", x$cluster, "s change an estimate by more than 1 SE;",
        " inspect them with plot().\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Tidy a leave-one-cluster-out diagnostic
#' @param x An object of class `lmmr_loso`.
#' @param ... Not used.
#' @return A data frame with one row per left-out cluster and coefficient.
#' @export
tidy.lmmr_loso <- function(x, ...) {
  data.frame(
    left_out = rep(rownames(x$estimates), ncol(x$estimates)),
    coef = rep(colnames(x$estimates), each = nrow(x$estimates)),
    estimate = as.vector(x$estimates),
    change = as.vector(x$change),
    std_change = as.vector(x$std_change),
    status = rep(x$status, ncol(x$estimates))
  )
}

#' Plot a leave-one-cluster-out diagnostic
#' @param x An object of class `lmmr_loso`.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object with one panel per coefficient.
#' @export
plot.lmmr_loso <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  td <- tidy.lmmr_loso(x)
  td$left_out <- factor(td$left_out, levels = rownames(x$estimates))
  ref <- data.frame(coef = names(x$full), full = unname(x$full),
                    lower = unname(x$full - 1.96 * x$se),
                    upper = unname(x$full + 1.96 * x$se))
  ggplot2::ggplot(td, ggplot2::aes(x = .data$estimate, y = .data$left_out)) +
    ggplot2::geom_rect(data = ref, inherit.aes = FALSE,
                       ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                                    ymin = -Inf, ymax = Inf),
                       fill = "grey90") +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = .data$full),
                        colour = "#C0392B") +
    ggplot2::geom_point(ggplot2::aes(colour = abs(.data$std_change) > 1),
                        size = 1.8) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20",
                                            `TRUE` = "#C0392B"),
                                 guide = "none") +
    ggplot2::facet_wrap(~coef, scales = "free_x") +
    ggplot2::labs(
      x = "Estimate without this cluster", y = paste0("Left-out ", x$cluster),
      title = "Leave-one-out influence",
      subtitle = "Red line: full-data estimate; band: its 95% Wald interval"
    ) +
    theme_lmmr()
}

#' Agreement between the model coefficient and participant-level differences
#'
#' Compares the model coefficient of a two-level within-cluster factor (e.g.
#' condition) with the mean of the differences computed directly in each
#' cluster (e.g. the mean of condition B minus the mean of condition A for
#' each participant). Both summarise the same effect, so a large discrepancy
#' suggests that the model is estimating something different from the
#' intended contrast.
#'
#' @details
#' A common cause in time-series models is a time course that differs between
#' conditions without a fixed effect of time or a condition by time
#' interaction in the model, so that the condition coefficient absorbs part
#' of the time course. Discrepancies can also arise legitimately: covariates
#' in the model, unequal numbers of trials per participant (the model weights
#' participants by their precision) or strong shrinkage of random slopes. The
#' check flags a discrepancy only when it exceeds both `tolerance` (relative
#' to the mean difference) and the standard error of the mean difference.
#'
#' @param model A linear model fitted with [stats::lm()] or [lme4::lmer()].
#' @param term Name of a two-level factor that varies within clusters.
#' @param cluster Name of the column identifying the clusters.
#' @param tolerance Relative discrepancy above which a warning is given.
#' @param data The data used to fit `model`, if it cannot be recovered from
#'   the model call.
#'
#' @return An object of class `lmmr_agreement` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()].
#' @seealso [diag_timecourse()]
#'
#' @examples
#' d <- sim_blocks("within", n_participants = 10, n_trials = 10, n_time = 20,
#'                 effect_condition = 0.5, seed = 1)
#' m <- lme4::lmer(y ~ condition + time + (1 | participant), data = d)
#' diag_agreement(m, term = "condition", cluster = "participant")
#' @export
diag_agreement <- function(model, term, cluster, tolerance = 0.25,
                           data = NULL) {
  check_model(model)
  check_column_name(term, "term")
  check_column_name(cluster, "cluster")
  check_number(tolerance, "tolerance", min = 0)
  data <- model_data(model, data)
  missing <- setdiff(c(term, cluster), names(data))
  if (length(missing) > 0) {
    cli::cli_abort("Column{?s} {.field {missing}} not found in the data.")
  }
  lev <- if (is.factor(data[[term]])) {
    levels(droplevels(data[[term]]))
  } else {
    sort(unique(as.character(data[[term]])))
  }
  if (length(lev) != 2) {
    cli::cli_abort("{.field {term}} must have exactly two levels.")
  }
  coef <- resolve_coef(model, term, NULL)
  response <- response_name(model)

  y <- data[[response]]
  grp <- as.character(data[[term]])
  cl <- as.character(data[[cluster]])
  means <- tapply(y, list(cl, grp), mean)
  diffs <- means[, lev[2]] - means[, lev[1]]
  diffs <- diffs[!is.na(diffs)]
  if (length(diffs) < 2) {
    cli::cli_abort("Fewer than two clusters have both levels of
                    {.field {term}}.")
  }
  mean_diff <- mean(diffs)
  se_diff <- stats::sd(diffs) / sqrt(length(diffs))
  estimate <- fixed_estimates(model)[[coef]]
  discrepancy <- estimate - mean_diff
  relative <- if (mean_diff != 0) discrepancy / abs(mean_diff) else Inf
  flagged <- abs(relative) > tolerance && abs(discrepancy) > se_diff

  structure(
    list(estimate = estimate, coef = coef, mean_diff = mean_diff,
         se_diff = se_diff, discrepancy = discrepancy, relative = relative,
         flagged = flagged, tolerance = tolerance, diffs = diffs,
         levels = lev, term = term, cluster = cluster),
    class = "lmmr_agreement"
  )
}

#' @export
print.lmmr_agreement <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  cat("\nAgreement for `", x$term, "` (", x$levels[2], " - ", x$levels[1],
      ")\n\n", sep = "")
  cat("Model coefficient:              ", f(x$estimate), "\n", sep = "")
  cat("Mean ", x$cluster, "-level difference: ", f(x$mean_diff),
      " (SE ", f(x$se_diff), ", ", length(x$diffs), " ", x$cluster, "s)\n",
      sep = "")
  cat("Discrepancy:                    ", f(x$discrepancy), " (",
      sprintf("%.0f", 100 * x$relative), "%)\n", sep = "")
  if (x$flagged) {
    cat("\nThe model coefficient differs from the ", x$cluster,
        "-level differences.\nCheck that the model includes the time course",
        " (fixed effect of time, and a\ncondition by time interaction if the",
        " curves differ); covariates or very\nunequal numbers of trials can",
        " also explain the difference.\n", sep = "")
  } else {
    cat("\nThe model coefficient agrees with the ", x$cluster,
        "-level differences.\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Tidy an agreement diagnostic
#' @param x An object of class `lmmr_agreement`.
#' @param ... Not used.
#' @return A data frame with one row per cluster and its difference.
#' @export
tidy.lmmr_agreement <- function(x, ...) {
  data.frame(cluster = names(x$diffs), difference = unname(x$diffs))
}

#' Plot an agreement diagnostic
#' @param x An object of class `lmmr_agreement`.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_agreement <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  td <- tidy.lmmr_agreement(x)
  td$cluster <- factor(td$cluster, levels = td$cluster[order(td$difference)])
  lines <- data.frame(
    value = c(x$mean_diff, x$estimate),
    what = c(paste0("Mean ", x$cluster, " difference"), "Model coefficient")
  )
  ggplot2::ggplot(td, ggplot2::aes(x = .data$difference, y = .data$cluster)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_point(colour = "grey25", size = 1.8) +
    ggplot2::geom_vline(data = lines,
                        ggplot2::aes(xintercept = .data$value,
                                     colour = .data$what,
                                     linetype = .data$what),
                        linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c("grey15", "#C0392B"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
    ggplot2::labs(
      x = paste0(x$levels[2], " - ", x$levels[1]), y = NULL,
      title = paste0("Agreement: ", x$term),
      subtitle = paste0("Differences computed within each ", x$cluster)
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
