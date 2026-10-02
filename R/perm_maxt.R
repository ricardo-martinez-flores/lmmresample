#' Specify one test for a family of permutation tests
#'
#' Pairs a fitted model with the term tested in it, for use in [perm_maxt()].
#'
#' @inheritParams perm_test
#' @return An object of class `lmmr_spec`.
#' @seealso [perm_maxt()]
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 6,
#'                 n_trials = 6, seed = 1)
#' perm_spec(lm(y ~ condition + participant, data = d), "condition")
#' @export
perm_spec <- function(model, term, coef = NULL, data = NULL) {
  check_model(model)
  check_column_name(term, "term")
  if (!is.null(coef) && (!is.character(coef) || anyNA(coef))) {
    cli::cli_abort("{.arg coef} must be a character vector of coefficient
                    names.")
  }
  data <- model_data(model, data)
  term_info(model, term)
  structure(list(model = model, term = term, coef = coef, data = data),
            class = "lmmr_spec")
}

#' @export
print.lmmr_spec <- function(x, ...) {
  cat("<test specification> `", x$term, "` in a ", class(x$model)[1],
      " model of `", deparse(stats::formula(x$model)[[2]]), "`\n", sep = "")
  invisible(x)
}

#' Family-wise permutation tests (max-t and min-p)
#'
#' Tests several model terms at once, for example the same condition effect
#' in one model per feature extracted from the same trials, and controls the
#' family-wise error rate with the max-t or min-p procedures of Westfall and
#' Young (1993). In each permutation the same resampling of units is applied to
#' every model, so the dependence among the test statistics is preserved and
#' the correction is less conservative than Bonferroni when the statistics
#' are correlated.
#'
#' @details
#' All models must be fitted to data that share the resampled units (e.g.
#' the same trials of the same participants). Two resampling methods are
#' available, as in [perm_test()]:
#'
#' * `"relabel"`: the values of the tested variable are reassigned among
#'   units. All specifications must then test the same variable, with the
#'   same values in every data set.
#' * `"freedman-lane"`: the residuals of each reduced model are permuted or
#'   sign-flipped with the same draw for every model. Specifications may test
#'   different terms, including interactions and effects adjusted for
#'   covariates.
#'
#' With `method = "auto"`, relabelling is used when it is valid for every
#' specification, and Freedman-Lane otherwise.
#'
#' Two ways of combining the tests are available:
#'
#' * `combine = "max-t"`: in each permutation the largest statistic across
#'   tests is recorded. The statistics must be on the same scale, so every
#'   term must have a single coefficient (absolute *t* statistics).
#' * `combine = "min-p"`: each statistic is first converted into a p-value
#'   within its own permutation distribution, and the smallest p-value across
#'   tests is recorded. It accepts statistics on different scales, such as a
#'   *t* for a main effect and a chi-square for a condition by spline-of-time
#'   interaction, and gives each test the same weight.
#'
#' With `combine = "auto"`, max-t is used when every term has a single
#' coefficient and min-p otherwise.
#'
#' * `adjust = "single-step"`: the adjusted p-value of each test is the
#'   proportion of permutations whose maximum statistic (or minimum p-value)
#'   is at least as extreme as the observed one for that test.
#' * `adjust = "step-down"`: tests are ordered from the most to the least
#'   extreme and the maximum (or minimum) is taken only over the tests not yet
#'   rejected, with monotonicity enforced. It is uniformly more powerful than
#'   the single-step procedure and controls the family-wise error rate under
#'   the same conditions.
#'
#' Adjusted p-values use the Phipson and Smyth (2010) correction. Permutations
#' in which any model failed to fit are excluded and counted. Tests are
#' two-sided.
#'
#' @param ... Named test specifications created with [perm_spec()]. Names are
#'   used to label the tests.
#' @param unit,exchange,time,B,seed As in [perm_test()].
#' @param method Resampling method: `"auto"`, `"relabel"` or
#'   `"freedman-lane"`. See Details.
#' @param combine `"auto"`, `"max-t"` or `"min-p"`. See Details.
#' @param adjust `"step-down"` (default) or `"single-step"`. See Details.
#'
#' @return An object of class `lmmr_maxt` with methods for [print()],
#'   [tidy()][generics::tidy] and [plot()]. Its element `results` is a data
#'   frame with the observed statistic, the unadjusted permutation p-value
#'   and the family-wise adjusted p-value of each test; `null` is the matrix
#'   of permuted statistics (one column per test).
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-based multiple
#' testing: Examples and methods for p-value adjustment*. Wiley.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero. *Statistical Applications in Genetics and Molecular Biology, 9*(1),
#' Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso [perm_spec()], [perm_test()]
#'
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 10,
#'                 n_trials = 8, effect_condition = 0.5, seed = 1)
#' m_mean <- lm(y ~ condition + participant, data = d)
#' m_peak <- lm(peak ~ condition + participant, data = d)
#'
#' # B is small here to keep the example fast; use the default for reporting
#' res <- perm_maxt(
#'   mean = perm_spec(m_mean, "condition"),
#'   peak = perm_spec(m_peak, "condition"),
#'   unit = "trial", exchange = exch_signflip("participant"),
#'   B = 99, seed = 1
#' )
#' res
#' @export
perm_maxt <- function(...,
                      unit,
                      exchange,
                      method = c("auto", "relabel", "freedman-lane"),
                      time = NULL,
                      combine = c("auto", "max-t", "min-p"),
                      adjust = c("step-down", "single-step"),
                      B = 4999,
                      seed = NULL) {
  call <- match.call()
  method <- match.arg(method)
  combine <- match.arg(combine)
  adjust <- match.arg(adjust)
  check_count(B, "B", min = 1)
  specs <- list(...)
  check_specs(specs)
  k <- length(specs)

  tests <- lapply(specs, function(s) {
    resolve_test(s$model, s$data, s$term, s$coef, method)
  })
  methods <- vapply(tests, function(t) t$method, character(1))
  if (method == "auto" && length(unique(methods)) > 1) {
    tests <- lapply(specs, function(s) {
      resolve_test(s$model, s$data, s$term, s$coef, "freedman-lane")
    })
  }
  resampling <- tests[[1]]$method
  types <- vapply(tests, function(t) t$type, character(1))
  if (combine == "auto") combine <- if (all(types == "t")) "max-t" else "min-p"
  if (combine == "max-t" && any(types != "t")) {
    cli::cli_abort(c(
      "With {.code combine = \"max-t\"} every tested term must have a single
       coefficient.",
      "i" = "Use {.code combine = \"min-p\"} to combine statistics on
             different scales."
    ))
  }
  if (resampling == "relabel") {
    terms <- unique(vapply(specs, function(s) s$term, character(1)))
    if (length(terms) > 1) {
      cli::cli_abort(c(
        "With relabelling, all specifications must test the same variable.",
        "x" = "Found {.val {terms}}.",
        "i" = "Use {.code method = \"freedman-lane\"} to test different
               terms."
      ))
    }
  }

  engines <- lapply(seq_len(k), function(j) {
    build_engine(specs[[j]]$model, specs[[j]]$data, specs[[j]]$term,
                 tests[[j]], unit, exchange, time)
  })
  ref <- engines[[1]]
  check_permutation_space(ref$space, B, exchange, unit)
  for (j in seq_len(k)[-1]) {
    check_same_engine_units(ref, engines[[j]], names(specs)[j], resampling)
  }

  observed <- vapply(seq_len(k), function(j) {
    if (types[j] == "chi2") {
      wald_chi2(specs[[j]]$model, tests[[j]]$coefs)[[1]]
    } else {
      abs(coef_stats(specs[[j]]$model, tests[[j]]$coefs)[[1]])
    }
  }, numeric(1))
  if (anyNA(observed)) {
    cli::cli_abort("The observed statistic is not available for
                    {.val {names(specs)[is.na(observed)]}}.")
  }

  draws_ref <- with_seed(seed, ref$draw(B))
  draws <- lapply(engines, function(e) e$translate(draws_ref, ref))
  refits <- lapply(specs, function(s) make_refitter(s$model))

  use_progress <- requireNamespace("progressr", quietly = TRUE)
  if (use_progress) p <- progressr::progressor(steps = B)
  results <- future.apply::future_lapply(seq_len(B), function(b) {
    out <- lapply(seq_len(k), function(j) {
      safe_refit(refits[[j]], engines[[j]]$make_data(draws[[j]], b),
                 tests[[j]]$coefs, types[j])
    })
    if (use_progress) p()
    out
  }, future.seed = FALSE,
  future.packages = do.call(refit_packages, lapply(specs, function(s) s$model)))

  null <- matrix(NA_real_, nrow = B, ncol = k,
                 dimnames = list(NULL, names(specs)))
  status <- matrix(NA_character_, nrow = B, ncol = k,
                   dimnames = list(NULL, names(specs)))
  for (b in seq_len(B)) {
    for (j in seq_len(k)) {
      r <- results[[b]][[j]]
      status[b, j] <- r$status
      if (!r$status %in% c("failed", "nonconverged")) {
        null[b, j] <- abs(r$stat[[1]])
      }
    }
  }
  complete <- stats::complete.cases(null)
  if (sum(!complete) / B > 0.05) {
    cli::cli_warn(c(
      "{sum(!complete)} of {B} permutations ({round(100 * sum(!complete) / B)}%)
       had at least one failed or non-convergent refit and were excluded.",
      "i" = "Consider simplifying the random-effects structure."
    ))
  }
  null_used <- null[complete, , drop = FALSE]

  p_unadj <- vapply(seq_len(k), function(j) {
    perm_pvalue(observed[j], null_used[, j], "greater")
  }, numeric(1))
  p_adj <- if (combine == "max-t") {
    maxt_adjust(observed, null_used, adjust)
  } else {
    minp_adjust(observed, null_used, adjust)
  }

  structure(
    list(
      results = data.frame(
        test = names(specs),
        term = vapply(specs, function(s) s$term, character(1)),
        coef = vapply(tests, function(t) paste(t$coefs, collapse = ", "),
                      character(1)),
        stat = ifelse(types == "chi2",
                      paste0("chi2(", vapply(tests, function(t) {
                        length(t$coefs)
                      }, integer(1)), ")"), "|t|"),
        statistic = unname(observed),
        p.value = p_unadj,
        p.adjusted = p_adj,
        row.names = NULL
      ),
      null = null,
      status = status,
      term = paste(unique(vapply(specs, function(s) s$term, character(1))),
                   collapse = ", "),
      unit = unit,
      exchange = exchange,
      method = resampling,
      adjust = adjust,
      combine = combine,
      stat_label = "t",
      n_units = ref$n_units,
      n_blocks = ref$space$n_blocks,
      B = B,
      B_used = sum(complete),
      seed = seed,
      call = call
    ),
    class = "lmmr_maxt"
  )
}

# All engines in a family must resample the same units with the same blocks.
check_same_engine_units <- function(ref, engine, name, resampling,
                                    call = rlang::caller_env()) {
  if (!setequal(ref$unit_ids, engine$unit_ids)) {
    cli::cli_abort(c(
      "Specification {.val {name}} does not have the same units as the
       first one.",
      "i" = "All models must be fitted to the same trials or participants."
    ), call = call)
  }
  idx <- match(ref$unit_ids, engine$unit_ids)
  if (!identical(ref$unit_blocks, engine$unit_blocks[idx])) {
    cli::cli_abort("Specification {.val {name}} assigns different blocks to
                    the same units.", call = call)
  }
  if (resampling == "relabel" &&
      !identical(ref$unit_values, engine$unit_values[idx])) {
    cli::cli_abort("Specification {.val {name}} assigns different values of
                    the tested variable to the same units.", call = call)
  }
  invisible(NULL)
}

# Westfall-Young max-t adjusted p-values (two-sided, absolute statistics).
maxt_adjust <- function(observed, null, method) {
  k <- length(observed)
  n <- nrow(null)
  if (n == 0) return(rep(NA_real_, k))
  tol <- 1e-8 * pmax(1, observed)
  if (method == "single-step") {
    max_null <- apply(null, 1, max)
    return(vapply(seq_len(k), function(j) {
      (sum(max_null >= observed[j] - tol[j]) + 1) / (n + 1)
    }, numeric(1)))
  }
  ord <- order(observed, decreasing = TRUE)
  p <- numeric(k)
  for (r in seq_len(k)) {
    j <- ord[r]
    remaining <- ord[r:k]
    max_null <- apply(null[, remaining, drop = FALSE], 1, max)
    p[j] <- (sum(max_null >= observed[j] - tol[j]) + 1) / (n + 1)
  }
  # Enforce monotonicity along the ordering
  p[ord] <- cummax(p[ord])
  p
}

# Permutation p-values of the observed statistics (row 1) and of every
# permuted statistic, each computed within its own column of
# rbind(observed, null), so that the observed data are part of the reference
# set (larger statistics are more extreme).
column_pvalues <- function(observed, null) {
  all <- rbind(observed, null)
  n <- nrow(all)
  matrix(apply(all, 2, function(v) (n - rank(v, ties.method = "min") + 1) / n),
         nrow = n)
}

# Westfall-Young min-p adjusted p-values.
minp_adjust <- function(observed, null, method) {
  k <- length(observed)
  if (nrow(null) == 0) return(rep(NA_real_, k))
  pv <- column_pvalues(observed, null)
  p_obs <- pv[1, ]
  tol <- 1e-12
  if (method == "single-step") {
    min_all <- apply(pv, 1, min)
    return(vapply(seq_len(k), function(j) {
      mean(min_all <= p_obs[j] + tol)
    }, numeric(1)))
  }
  ord <- order(p_obs)
  p <- numeric(k)
  for (r in seq_len(k)) {
    j <- ord[r]
    remaining <- ord[r:k]
    min_all <- apply(pv[, remaining, drop = FALSE], 1, min)
    p[j] <- mean(min_all <= p_obs[j] + tol)
  }
  p[ord] <- cummax(p[ord])
  p
}

check_specs <- function(specs, call = rlang::caller_env()) {
  if (length(specs) < 2) {
    cli::cli_abort("{.fn perm_maxt} needs at least two specifications.",
                   call = call)
  }
  if (!all(vapply(specs, inherits, logical(1), "lmmr_spec"))) {
    cli::cli_abort("All arguments in {.arg ...} must be created with
                    {.fn perm_spec}.", call = call)
  }
  nms <- names(specs)
  if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    cli::cli_abort("Specifications must have unique names, e.g.
                    {.code perm_maxt(peak = perm_spec(...), latency =
                    perm_spec(...), ...)}.", call = call)
  }
  invisible(specs)
}

# Methods -------------------------------------------------------------------

#' @export
print.lmmr_maxt <- function(x, digits = 3, ...) {
  cat("\nFamily-wise permutation tests for `", x$term, "` (", x$combine,
      ", ", x$adjust, ")\n\n", sep = "")
  cat("Method:    ", method_label(x$method), "\n", sep = "")
  cat("Exchange:  ", format(x$exchange, unit = x$unit), " (", x$n_units,
      " units", if (!is.null(x$exchange$block)) {
        paste0(" in ", x$n_blocks, " blocks")
      }, ")\n", sep = "")
  cat("Resamples: ", x$B_used, " of ", x$B, " permutations used\n\n",
      sep = "")
  res <- x$results
  tab <- data.frame(
    Test = res$test,
    Statistic = paste(res$stat, "=",
                      formatC(res$statistic, digits = digits, format = "f")),
    `p (unadjusted)` = vapply(res$p.value, format_p, character(1)),
    `p (family-wise)` = vapply(res$p.adjusted, format_p, character(1)),
    check.names = FALSE
  )
  print(tab, row.names = FALSE, right = FALSE)
  if (x$B < 1000) {
    cat("\nNote: fewer than 1000 permutations; treat these results as",
        "exploratory.\n")
  }
  cat("\n")
  invisible(x)
}

#' Tidy family-wise permutation tests
#'
#' @param x An object of class `lmmr_maxt`.
#' @param type `"summary"` returns one row per test; `"null"` returns the
#'   permuted statistics in long format.
#' @param ... Not used.
#' @return A data frame.
#' @export
tidy.lmmr_maxt <- function(x, type = c("summary", "null"), ...) {
  type <- match.arg(type)
  if (type == "null") {
    return(data.frame(
      permutation = rep(seq_len(nrow(x$null)), ncol(x$null)),
      test = rep(colnames(x$null), each = nrow(x$null)),
      statistic = as.vector(x$null)
    ))
  }
  out <- x$results
  out$method <- paste0(x$combine, " (", x$adjust, "), ", x$method)
  out
}

#' Plot family-wise permutation tests
#'
#' Shows the permutation distribution of the maximum absolute statistic
#' (max-t) or of the minimum p-value (min-p) across tests, with the observed
#' value of each test and the critical value that controls the family-wise
#' error rate at `alpha`.
#'
#' @param x An object of class `lmmr_maxt`.
#' @param alpha Family-wise significance level for the critical value.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_maxt <- function(x, alpha = 0.05, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  null <- x$null[stats::complete.cases(x$null), , drop = FALSE]
  obs <- x$results
  obs$significant <- obs$p.adjusted <= alpha
  if (identical(x$combine, "min-p")) {
    pv <- column_pvalues(obs$statistic, null)
    min_p <- apply(pv[-1, , drop = FALSE], 1, min)
    critical <- stats::quantile(min_p, alpha, names = FALSE)
    return(
      ggplot2::ggplot(data.frame(value = min_p),
                      ggplot2::aes(x = .data$value)) +
        ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = "white",
                                linewidth = 0.2) +
        ggplot2::geom_vline(xintercept = critical, linetype = "dashed",
                            colour = "grey30") +
        ggplot2::geom_vline(
          data = obs,
          ggplot2::aes(xintercept = .data$p.value, colour = .data$significant),
          linewidth = 0.8
        ) +
        ggplot2::geom_text(
          data = obs,
          ggplot2::aes(x = .data$p.value, y = Inf, label = .data$test,
                       colour = .data$significant),
          angle = 90, hjust = 1.1, vjust = -0.4, size = 3.2,
          show.legend = FALSE
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(
          values = c(`TRUE` = "#C0392B", `FALSE` = "grey20"),
          labels = c(`TRUE` = "significant", `FALSE` = "not significant"),
          name = NULL
        ) +
        ggplot2::labs(
          x = "Minimum p-value across tests (log scale)", y = "Permutations",
          title = "Family-wise null distribution (min-p)",
          subtitle = paste0("Lines: unadjusted p of each test; dashed: ",
                            "critical value at family-wise alpha = ", alpha)
        ) +
        theme_lmmr() +
        ggplot2::theme(legend.position = "bottom")
    )
  }
  max_null <- apply(null, 1, max)
  critical <- stats::quantile(max_null, 1 - alpha, names = FALSE)
  ggplot2::ggplot(data.frame(value = max_null),
                  ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 50, fill = "grey75", colour = "white",
                            linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = critical, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::geom_vline(
      data = obs,
      ggplot2::aes(xintercept = .data$statistic, colour = .data$significant),
      linewidth = 0.8
    ) +
    ggplot2::geom_text(
      data = obs,
      ggplot2::aes(x = .data$statistic, y = Inf, label = .data$test,
                   colour = .data$significant),
      angle = 90, hjust = 1.1, vjust = -0.4, size = 3.2, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#C0392B", `FALSE` = "grey20"),
      labels = c(`TRUE` = "significant", `FALSE` = "not significant"),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Maximum |t| across tests",
      y = "Permutations",
      title = paste0("Family-wise null distribution: ", x$term),
      subtitle = paste0("Dashed line: critical value at family-wise alpha = ",
                        alpha, " (", round(critical, 2), "); ", x$B_used,
                        " permutations")
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
