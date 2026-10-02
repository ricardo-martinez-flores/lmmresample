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
  if (!is.null(coef)) check_column_name(coef, "coef")
  data <- model_data(model, data)
  coef <- resolve_coef(model, term, coef)
  structure(list(model = model, term = term, coef = coef, data = data),
            class = "lmmr_spec")
}

#' @export
print.lmmr_spec <- function(x, ...) {
  cat("<test specification> `", x$coef, "` in a ", class(x$model)[1],
      " model of `", deparse(stats::formula(x$model)[[2]]), "`\n", sep = "")
  invisible(x)
}

#' Family-wise permutation tests with the max-t procedure
#'
#' Tests the same variable in several models at once, for example one model
#' per feature extracted from the same trials, and controls the family-wise
#' error rate with the max-t procedure of Westfall and Young (1993). In each
#' permutation the same relabelling of units is applied to every model, so
#' the dependence among the test statistics is preserved and the correction
#' is less conservative than Bonferroni when the statistics are correlated.
#'
#' @details
#' All specifications must test the same variable (`term`) in data sets that
#' share the units and their values of `term`, as when several outcome
#' variables are measured on the same trials. For each permutation, the
#' absolute Wald statistics of all models are computed and their maximum is
#' recorded.
#'
#' * `method = "single-step"`: the adjusted p-value of each test is the
#'   proportion of permutations whose maximum statistic is at least as large
#'   as the observed statistic of that test.
#' * `method = "step-down"`: tests are ordered by decreasing observed
#'   statistic and the maximum is taken only over the tests not yet rejected,
#'   with monotonicity enforced. It is uniformly more powerful than the
#'   single-step procedure and controls the family-wise error rate under the
#'   same conditions.
#'
#' Adjusted p-values use the Phipson and Smyth (2010) correction. Permutations
#' in which any model failed to fit are excluded and counted. Tests are
#' two-sided.
#'
#' @param ... Named test specifications created with [perm_spec()]. Names are
#'   used to label the tests.
#' @param unit,exchange,B,seed As in [perm_test()].
#' @param method `"step-down"` (default) or `"single-step"`. See Details.
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
                      method = c("step-down", "single-step"),
                      B = 4999,
                      seed = NULL) {
  call <- match.call()
  method <- match.arg(method)
  check_count(B, "B", min = 1)
  specs <- list(...)
  check_specs(specs)
  term <- specs[[1]]$term

  # Units of the first specification define the common relabelling
  ub <- lapply(specs, function(s) build_units(s$data, term, unit, exchange))
  ref <- ub[[1]]$units
  values <- specs[[1]]$data[[term]][ref$value_row]
  for (j in seq_along(specs)[-1]) {
    check_same_units(ref, values, ub[[j]]$units,
                     specs[[j]]$data[[term]][ub[[j]]$units$value_row],
                     names(specs)[j])
  }
  enc <- encode_values(values, exchange)
  n_perm <- count_permutations(values, ref$block, exchange$type)
  check_permutation_space(n_perm, B, exchange, unit)

  # Map each specification's rows to the common unit order
  row_map <- lapply(ub, function(u) match(u$units$unit, ref$unit)[u$row_unit])

  observed <- vapply(specs, function(s) abs(coef_stats(s$model, s$coef)[[1]]),
                     numeric(1))
  if (anyNA(observed)) {
    cli::cli_abort("The observed statistic is not available for
                    {.val {names(specs)[is.na(observed)]}}.")
  }

  relabels <- with_seed(seed, generate_relabels(enc$codes, ref$block,
                                                 exchange$type, B))
  distinct <- enc$distinct
  refits <- lapply(specs, function(s) make_refitter(s$model))

  use_progress <- requireNamespace("progressr", quietly = TRUE)
  if (use_progress) p <- progressr::progressor(steps = B)
  results <- future.apply::future_lapply(seq_len(B), function(b) {
    new_values <- distinct[relabels[, b]]
    out <- lapply(seq_along(specs), function(j) {
      newdata <- specs[[j]]$data
      newdata[[term]] <- new_values[row_map[[j]]]
      safe_refit(refits[[j]], newdata, specs[[j]]$coef)
    })
    if (use_progress) p()
    out
  }, future.seed = FALSE)

  k <- length(specs)
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
    perm_pvalue(observed[j], null_used[, j], "two.sided")
  }, numeric(1))
  p_adj <- maxt_adjust(observed, null_used, method)

  structure(
    list(
      results = data.frame(
        test = names(specs),
        coef = vapply(specs, function(s) s$coef, character(1)),
        statistic = unname(observed),
        p.value = p_unadj,
        p.adjusted = p_adj,
        row.names = NULL
      ),
      null = null,
      status = status,
      term = term,
      unit = unit,
      exchange = exchange,
      method = method,
      stat_label = stat_label(specs[[1]]$model),
      n_units = nrow(ref),
      n_blocks = n_perm$n_blocks,
      B = B,
      B_used = sum(complete),
      seed = seed,
      call = call
    ),
    class = "lmmr_maxt"
  )
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
  terms <- unique(vapply(specs, function(s) s$term, character(1)))
  if (length(terms) > 1) {
    cli::cli_abort(c(
      "All specifications must test the same variable.",
      "x" = "Found {.val {terms}}."
    ), call = call)
  }
  invisible(specs)
}

check_same_units <- function(ref, values, units, unit_values, name,
                             call = rlang::caller_env()) {
  if (!setequal(ref$unit, units$unit)) {
    cli::cli_abort(c(
      "Specification {.val {name}} does not have the same units as the
       first one.",
      "i" = "All models must be fitted to the same trials or participants."
    ), call = call)
  }
  idx <- match(ref$unit, units$unit)
  if (!identical(as.character(values), as.character(unit_values[idx])) ||
      !identical(ref$block, units$block[idx])) {
    cli::cli_abort(
      "Specification {.val {name}} assigns different values of the tested
       variable or different blocks to the same units.",
      call = call
    )
  }
  invisible(NULL)
}

# Methods -------------------------------------------------------------------

#' @export
print.lmmr_maxt <- function(x, digits = 3, ...) {
  cat("\nFamily-wise permutation tests for `", x$term, "` (max-t, ",
      x$method, ")\n\n", sep = "")
  cat("Exchange:  ", format(x$exchange, unit = x$unit), " (", x$n_units,
      " units", if (!is.null(x$exchange$block)) {
        paste0(" in ", x$n_blocks, " blocks")
      }, ")\n", sep = "")
  cat("Resamples: ", x$B_used, " of ", x$B, " permutations used\n\n",
      sep = "")
  res <- x$results
  tab <- data.frame(
    Test = res$test,
    Statistic = formatC(res$statistic, digits = digits, format = "f"),
    `p (unadjusted)` = vapply(res$p.value, format_p, character(1)),
    `p (family-wise)` = vapply(res$p.adjusted, format_p, character(1)),
    check.names = FALSE
  )
  names(tab)[2] <- paste0("|", x$stat_label, "|")
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
  out$method <- paste0("max-t (", x$method, ")")
  out
}

#' Plot family-wise permutation tests
#'
#' Shows the permutation distribution of the maximum absolute statistic
#' across tests, with the observed statistic of each test and the critical
#' value that controls the family-wise error rate at `alpha`.
#'
#' @param x An object of class `lmmr_maxt`.
#' @param alpha Family-wise significance level for the critical value.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.lmmr_maxt <- function(x, alpha = 0.05, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  null <- x$null[stats::complete.cases(x$null), , drop = FALSE]
  max_null <- apply(null, 1, max)
  critical <- stats::quantile(max_null, 1 - alpha, names = FALSE)
  obs <- x$results
  obs$significant <- obs$p.adjusted <= alpha
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
      x = paste0("Maximum |", x$stat_label, "| across tests"),
      y = "Permutations",
      title = paste0("Family-wise null distribution: ", x$term),
      subtitle = paste0("Dashed line: critical value at family-wise alpha = ",
                        alpha, " (", round(critical, 2), "); ", x$B_used,
                        " permutations")
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
