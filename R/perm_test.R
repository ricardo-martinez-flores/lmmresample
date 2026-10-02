#' Permutation test for a model term
#'
#' Tests a fixed-effect term by resampling whole units (trials or
#' participants) and refitting the full model to each resampled data set.
#' Because complete units are exchanged, any dependence within a unit, such
#' as the autocorrelation of samples within a trial, is preserved in the
#' permutation distribution without having to be modelled.
#'
#' @details
#' Two resampling methods are available:
#'
#' * `"relabel"`: the values of `term` are reassigned among units, as
#'   declared by `exchange`. It requires `term` to be a variable constant
#'   within units and is used for main effects.
#' * `"freedman-lane"`: the reduced model without `term` is fitted, and its
#'   residuals are permuted between whole units (or sign-flipped by block with
#'   [exch_signflip()]) and added back to its fitted values (Freedman & Lane,
#'   1983; Winkler et al., 2014). It accounts for nuisance terms and allows
#'   tests of interactions (e.g. `"condition:time"`) and of terms with several
#'   coefficients. Random slopes of `term` are removed from the reduced model,
#'   because they would absorb the tested effect. Permuting residuals between units requires units of equal
#'   size within blocks; sign-flipping accepts units of any size.
#'
#' With `method = "auto"` (default), relabelling is used for a main effect that
#' is a column of the data with a single coefficient, and Freedman-Lane
#' otherwise.
#'
#' The statistic is the Wald statistic of the coefficient (its *t* value:
#' estimate divided by its standard error). When the term has several coefficients
#' (a factor with more than two levels, or an interaction with a spline of
#' time), the statistic is the Wald chi-square of all of them and the test is
#' two-sided. Refits keep every argument of the original call; only the data
#' change.
#'
#' The p-value follows Phipson and Smyth (2010):
#' \deqn{p = (b + 1) / (B' + 1),}
#' where \eqn{b} is the number of permutations with a statistic at least as
#' extreme as the observed one and \eqn{B'} is the number of permutations
#' that were refitted successfully. Refits that fail or do not converge are
#' excluded and counted; singular fits are kept and counted.
#'
#' Permutations are generated in the main R session, so results are
#' identical for a given `seed` whatever the parallel backend. Refits run in
#' parallel through the \pkg{future} framework: choose a backend with
#' [future::plan()], e.g. `future::plan("multisession", workers = 4)`.
#' Progress can be reported with [progressr::with_progress()].
#'
#' Terms contained in a higher-order interaction are refused: test the
#' interaction first.
#'
#' @param model A linear model fitted with [stats::lm()] or [lme4::lmer()].
#' @param term Fixed-effect term to test, as labelled in the model formula
#'   (e.g. `"condition"` or `"condition:time"`). For `method = "relabel"` it
#'   must be a column of the data, constant within each `unit`.
#' @param unit Name of the column identifying the units whose values of
#'   `term` are permuted (e.g. `"trial"` or `"participant"`).
#' @param exchange Exchangeability of units, created with [exch_within()],
#'   [exch_signflip()] or [exch_free()].
#' @param coef Optional name of the coefficient(s) of `term` to test. By
#'   default all coefficients of the term are tested jointly.
#' @param method Resampling method: `"auto"`, `"relabel"` or
#'   `"freedman-lane"`. See Details.
#' @param time Name of the column giving the order of samples within units.
#'   Required when Freedman-Lane residuals are permuted between units of
#'   several rows (with [exch_within()] or [exch_free()]), so that samples are
#'   matched by time; not needed with [exch_signflip()].
#' @param alternative Direction of the alternative hypothesis.
#' @param B Number of permutations. The default, 4999, is suitable for
#'   reporting; use smaller values only for exploration.
#' @param seed Optional integer seed for reproducible permutations. The global
#'   random number generator state is left unchanged.
#' @param data The data used to fit `model`. Needed only if it cannot be
#'   recovered from the model call.
#'
#' @return An object of class `lmmr_perm` with methods for [print()],
#'   [summary()], [tidy()][generics::tidy] and [plot()]. Its elements include:
#'   \item{statistic}{Observed statistic.}
#'   \item{null}{Statistics of all `B` permutations (`NA` for excluded
#'     refits).}
#'   \item{status}{Refit status of each permutation: `"ok"`, `"singular"`,
#'     `"nonconverged"` or `"failed"`.}
#'   \item{p.value}{Permutation p-value.}
#'   \item{B, B_used}{Requested and successfully used permutations.}
#'
#' @references
#' Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of
#' reported significance levels. *Journal of Business & Economic Statistics,
#' 1*(4), 292-298. \doi{10.1080/07350015.1983.10509354}
#'
#' Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols,
#' T. E. (2014). Permutation inference for the general linear model.
#' *NeuroImage, 92*, 381-397. \doi{10.1016/j.neuroimage.2014.01.060}
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: Calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology, 9*(1),
#' Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso [exch_within()], [exch_signflip()], [exch_free()],
#'   [perm_calibrate()]
#'
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 12,
#'                 n_trials = 10, effect_condition = 0.4, seed = 1)
#' m <- lm(y ~ condition + participant, data = d)
#'
#' # B is small here to keep the example fast; use the default for reporting
#' res <- perm_test(m, term = "condition", unit = "trial",
#'                  exchange = exch_within("participant"), B = 199, seed = 1)
#' res
#'
#' \donttest{
#' # Mixed model on the full time series
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   ts <- sim_blocks("within", n_participants = 10, n_trials = 10,
#'                    n_time = 50, effect_condition = 0.4, seed = 2)
#'   m_ts <- lme4::lmer(y ~ condition + (1 | participant), data = ts)
#'   perm_test(m_ts, term = "condition", unit = "trial",
#'             exchange = exch_within("participant"), B = 99, seed = 1)
#'
#'   # When does the effect appear? Condition x time interaction
#'   ts2 <- sim_blocks("within", n_participants = 10, n_trials = 10,
#'                     n_time = 30, effect_condition = 0.5,
#'                     effect_onset = 0.3, seed = 3)
#'   m_int <- lme4::lmer(y ~ condition * time + (1 | participant), data = ts2)
#'   perm_test(m_int, term = "condition:time", unit = "trial",
#'             exchange = exch_signflip("participant"), B = 99, seed = 1)
#' }
#' }
#' @export
perm_test <- function(model,
                      term,
                      unit,
                      exchange,
                      coef = NULL,
                      method = c("auto", "relabel", "freedman-lane"),
                      time = NULL,
                      alternative = c("two.sided", "greater", "less"),
                      B = 4999,
                      seed = NULL,
                      data = NULL) {
  call <- match.call()
  check_model(model)
  check_column_name(term, "term")
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  check_count(B, "B", min = 1)
  if (!is.null(coef) && (!is.character(coef) || anyNA(coef))) {
    cli::cli_abort("{.arg coef} must be a character vector of coefficient
                    names.")
  }

  data <- model_data(model, data)
  test <- resolve_test(model, data, term, coef, method)
  if (test$type == "chi2" && alternative != "two.sided") {
    cli::cli_abort("One-sided alternatives are only available for terms with
                    a single coefficient.")
  }
  engine <- build_engine(model, data, term, test, unit, exchange, time)
  check_permutation_space(engine$space, B, exchange, unit)

  observed <- if (test$type == "chi2") {
    wald_chi2(model, test$coefs)[[1]]
  } else {
    coef_stats(model, test$coefs)[[1]]
  }
  if (is.na(observed)) {
    cli::cli_abort("The observed statistic for {.val {term}} is not
                    available (coefficient not estimable).")
  }

  draws <- with_seed(seed, engine$draw(B))
  refit <- make_refitter(model)
  res <- run_refits(refit, function(b) engine$make_data(draws, b), B,
                    test$coefs, test$type)

  null <- res$stats[, 1]
  status <- res$status
  null[status %in% c("failed", "nonconverged")] <- NA_real_
  used <- !is.na(null)
  if (sum(!used) / B > 0.05) {
    cli::cli_warn(c(
      "{sum(!used)} of {B} refits ({round(100 * sum(!used) / B)}%) failed or
       did not converge and were excluded.",
      "i" = "The permutation distribution may be biased. Consider
             simplifying the random-effects structure."
    ))
  }
  p_alt <- if (test$type == "chi2") "greater" else alternative

  structure(
    list(
      statistic = observed,
      null = null,
      status = status,
      p.value = perm_pvalue(observed, null[used], p_alt),
      alternative = alternative,
      term = term,
      coef = test$coefs,
      method = test$method,
      stat_type = test$type,
      stat_label = if (test$type == "chi2") "chi2" else stat_label(model),
      df = length(test$coefs),
      unit = unit,
      exchange = exchange,
      n_units = engine$n_units,
      n_blocks = engine$space$n_blocks,
      n_uninformative = engine$space$uninformative,
      log_n_permutations = engine$space$log_n,
      B = B,
      B_used = sum(used),
      seed = seed,
      model_class = class(model)[1],
      call = call
    ),
    class = "lmmr_perm"
  )
}

perm_pvalue <- function(observed, null, alternative) {
  if (length(null) == 0) return(NA_real_)
  tol <- 1e-8 * max(1, abs(observed))
  b <- switch(alternative,
    two.sided = sum(abs(null) >= abs(observed) - tol),
    greater = sum(null >= observed - tol),
    less = sum(null <= observed + tol)
  )
  (b + 1) / (length(null) + 1)
}

check_permutation_space <- function(n_perm, B, exchange, unit,
                                    call = rlang::caller_env()) {
  if (n_perm$uninformative == n_perm$n_blocks) {
    cli::cli_abort(c(
      "No permutation can change the data.",
      "i" = "The tested variable takes a single value within every
             {if (exchange$type == 'free') 'stratum' else 'block'}."
    ), call = call)
  }
  if (n_perm$uninformative > 0 && exchange$type == "signflip") {
    cli::cli_abort(c(
      "{n_perm$uninformative} of {n_perm$n_blocks} block{?s} contain{?s/} a
       single value of the tested variable.",
      "i" = "{.fn exch_signflip} requires both values in every block."
    ), call = call)
  }
  if (n_perm$uninformative > 0 && exchange$type == "within") {
    cli::cli_warn(c(
      "{n_perm$uninformative} of {n_perm$n_blocks} block{?s} contain{?s/} a
       single value of the tested variable and do not contribute to the
       permutation distribution."
    ), call = call)
  }
  if (n_perm$log_n < log(B)) {
    total <- round(exp(n_perm$log_n))
    cli::cli_warn(c(
      "Only {total} distinct relabelling{?s} of {.field {unit}} exist{?s/},
       fewer than {.arg B} = {B}.",
      "i" = "Permutations are drawn with replacement; the test is valid but
             can be conservative. The smallest attainable p-value is about
             {signif(1 / total, 2)}."
    ), call = call)
  }
  invisible(NULL)
}

# Methods -------------------------------------------------------------------

#' @export
print.lmmr_perm <- function(x, digits = 3, ...) {
  cat("\nPermutation test for `", x$term, "`", sep = "")
  if (length(x$coef) > 1) {
    cat(" (", length(x$coef), " coefficients)", sep = "")
  } else if (x$coef != x$term) {
    cat(" (coefficient `", x$coef, "`)", sep = "")
  }
  cat("\n\n")
  cat("Method:    ", method_label(x$method), "\n", sep = "")
  cat("Exchange:  ", format(x$exchange, unit = x$unit), " (",
      x$n_units, " units", if (!is.null(x$exchange$block)) {
        paste0(" in ", x$n_blocks, " blocks")
      }, ")\n", sep = "")
  cat("Statistic: ", x$stat_label, if (x$stat_type == "chi2") {
    paste0("(", x$df, ")")
  }, " = ", formatC(x$statistic, digits = digits, format = "f"), "\n",
  sep = "")
  cat("p-value:   ", format_p(x$p.value), " (",
      if (x$stat_type == "chi2") "joint test" else alt_label(x$alternative),
      ")\n", sep = "")
  cat("Resamples: ", x$B_used, " of ", x$B, " permutations used\n", sep = "")
  print_refit_notes(x$status, x$B, x$B_used)
  cat("\n")
  invisible(x)
}

#' @export
summary.lmmr_perm <- function(object, ...) {
  structure(object, class = c("summary.lmmr_perm", class(object)))
}

#' @export
print.summary.lmmr_perm <- function(x, digits = 3, ...) {
  print.lmmr_perm(x, digits = digits)
  null <- x$null[!is.na(x$null)]
  cat("Permutation distribution of ", x$stat_label, ":\n", sep = "")
  q <- stats::quantile(null, c(0, 0.025, 0.5, 0.975, 1))
  names(q) <- c("min", "2.5%", "median", "97.5%", "max")
  print(round(q, digits))
  cat("\nRefit status:\n")
  print(table(factor(x$status,
                     levels = c("ok", "singular", "nonconverged", "failed"))))
  cat("\nDistinct permutations available: ",
      format_count(x$log_n_permutations), "\n", sep = "")
  if (!is.null(x$seed)) cat("Seed: ", x$seed, "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Tidy a permutation test
#'
#' @param x An object of class `lmmr_perm`.
#' @param type `"summary"` returns one row with the test result; `"null"`
#'   returns the permutation distribution, one row per permutation.
#' @param ... Not used.
#' @return A data frame.
#' @export
tidy.lmmr_perm <- function(x, type = c("summary", "null"), ...) {
  type <- match.arg(type)
  if (type == "null") {
    return(data.frame(permutation = seq_along(x$null), statistic = x$null,
                      status = x$status))
  }
  data.frame(
    term = x$term, coef = paste(x$coef, collapse = ", "),
    statistic = x$statistic, df = x$df, p.value = x$p.value,
    alternative = if (x$stat_type == "chi2") "joint" else x$alternative,
    B = x$B, B_used = x$B_used,
    method = paste0("permutation (", x$method, ")")
  )
}

#' Plot a permutation test
#'
#' @param x An object of class `lmmr_perm`.
#' @param type `"null"` shows the permutation distribution with the observed
#'   statistic. `"trace"` shows the running p-value as permutations
#'   accumulate, to check that `B` is large enough for a stable result.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_perm <- function(x, type = c("null", "trace"), ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  type <- match.arg(type)
  null <- x$null[!is.na(x$null)]
  stat_name <- paste0(x$stat_label, " statistic")

  if (type == "null") {
    obs <- data.frame(value = x$statistic, kind = "observed")
    if (x$alternative == "two.sided" && x$stat_type != "chi2") {
      obs <- rbind(obs, data.frame(value = -x$statistic, kind = "mirror"))
    }
    return(
      ggplot2::ggplot(data.frame(value = null), ggplot2::aes(x = .data$value)) +
        ggplot2::geom_histogram(bins = 50, fill = "grey75", colour = "white",
                                linewidth = 0.2) +
        ggplot2::geom_vline(
          data = obs,
          ggplot2::aes(xintercept = .data$value, linetype = .data$kind),
          colour = "#C0392B", linewidth = 0.8, show.legend = FALSE
        ) +
        ggplot2::scale_linetype_manual(
          values = c(observed = "solid", mirror = "dashed")
        ) +
        ggplot2::labs(
          x = stat_name, y = "Permutations",
          title = paste0("Permutation distribution: ", x$term),
          subtitle = paste0(x$stat_label, " = ", round(x$statistic, 2), ", p ",
                            format_p(x$p.value, prefix = TRUE), ", ",
                            x$B_used, " permutations")
        ) +
        theme_lmmr()
    )
  }

  alt <- if (x$stat_type == "chi2") "greater" else x$alternative
  extreme <- switch(alt,
    two.sided = abs(null) >= abs(x$statistic),
    greater = null >= x$statistic,
    less = null <= x$statistic
  )
  k <- seq_along(null)
  trace <- data.frame(k = k, p = (cumsum(extreme) + 1) / (k + 1))
  se <- sqrt(x$p.value * (1 - x$p.value) / k)
  trace$lower <- pmax(0, x$p.value - 1.96 * se)
  trace$upper <- pmin(1, x$p.value + 1.96 * se)
  ggplot2::ggplot(trace, ggplot2::aes(x = .data$k, y = .data$p)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_line(colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = x$p.value, colour = "#C0392B",
                        linetype = "dashed") +
    ggplot2::labs(
      x = "Number of permutations", y = "Running p-value",
      title = paste0("Stability of the p-value: ", x$term),
      subtitle = "Band: 95% Monte Carlo interval around the final p-value"
    ) +
    theme_lmmr()
}

# Formatting helpers ----------------------------------------------------------

format_p <- function(p, prefix = FALSE) {
  if (is.na(p)) return(if (prefix) "= NA" else "NA")
  out <- if (p < 0.001) "< 0.001" else formatC(p, digits = 3, format = "f")
  if (prefix && p >= 0.001) paste("=", out) else out
}

method_label <- function(method) {
  switch(method, relabel = "relabelling of units",
         `freedman-lane` = "Freedman-Lane (residuals of the reduced model)")
}

alt_label <- function(alternative) {
  switch(alternative, two.sided = "two-sided", greater = "one-sided, greater",
         less = "one-sided, less")
}

format_count <- function(log_n) {
  if (log_n < log(1e6)) return(format(round(exp(log_n)), big.mark = ","))
  paste0("about 10^", floor(log_n / log(10)))
}

print_refit_notes <- function(status, B, B_used) {
  n_singular <- sum(status == "singular")
  n_excluded <- B - B_used
  if (n_excluded > 0) {
    cat("           ", n_excluded, " refit", if (n_excluded > 1) "s",
        " excluded (failed or not converged)\n", sep = "")
  }
  if (n_singular > 0) {
    cat("           ", n_singular, " singular fit", if (n_singular > 1) "s",
        " kept\n", sep = "")
  }
  if (B < 1000) {
    cat("Note: fewer than 1000 permutations; treat this result as",
        "exploratory.\n")
  }
}

theme_lmmr <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35")
    )
}
