#' Pairwise comparisons between the levels of a factor
#'
#' Compares pairs of levels of a factor (post hoc tests) with Freedman-Lane
#' permutation tests, controls the family-wise error rate across pairs, and
#' optionally adds bootstrap confidence intervals of each pairwise
#' difference.
#'
#' @details
#' For each pair, the model is refitted with the factor coded by indicator
#' (dummy) variables whose reference is the first level of the pair; this
#' reparametrisation leaves the fit unchanged. The pair is tested through
#' the indicator of its second level, while the indicators of the other
#' levels remain in the reduced model of the Freedman-Lane procedure. A
#' difference involving other levels therefore does not leak into the test
#' of the pair, as it would if all labels of the factor were permuted
#' together.
#'
#' `term` sets what is compared: the factor itself (difference in level) or
#' an interaction that contains it, such as `"condition:ns(time, 3)"`
#' (difference in time course). Random slopes of the factor are
#' reparametrised in the same way.
#'
#' All pairs are tested on the same data with the same permutations, so they
#' form one family, adjusted by default with the step-down procedure of
#' Westfall and Young (1993): max-t when each pair is a single coefficient,
#' min-p otherwise (see [perm_maxt()]). Holm's method on the unadjusted
#' permutation p-values is also available.
#'
#' Confidence intervals (when `cluster` is given) are bootstrap intervals of
#' each pairwise difference from [boot_ci()], with case resampling of
#' `cluster`; they are available when each pair is a single coefficient.
#' They are not adjusted for multiplicity.
#'
#' @param model A model fitted with [stats::lm()] or [lme4::lmer()] that
#'   contains `factor`.
#' @param factor Name of the factor (or character column) whose levels are
#'   compared.
#' @param term The model term compared between levels: `factor` (default) or
#'   an interaction containing it.
#' @param unit,exchange,time,B,seed As in [perm_test()].
#' @param pairs Pairs of levels to compare: a two-column character matrix
#'   (or a list of length-two character vectors). Each row gives the
#'   reference level and the compared level. Default: all pairs, in the order
#'   of the factor levels.
#' @param adjust Family-wise adjustment: `"westfall-young"` (default; max-t
#'   or min-p, step-down), `"holm"` or `"none"`.
#' @param cluster Optional cluster column for bootstrap confidence intervals
#'   (e.g. `"participant"`). `NULL` (default) skips the intervals.
#' @param strata,ci,level Passed to [boot_ci()].
#' @param data The data used to fit `model`, if it cannot be recovered from
#'   the model call.
#'
#' @return An object of class `lmmr_pairwise` with methods for [print()] and
#'   [tidy()][generics::tidy]. Its element `results` is a data frame with one
#'   row per pair: the levels, the estimated difference (compared minus
#'   reference level; single-coefficient comparisons only), the statistic,
#'   the unadjusted and adjusted permutation p-values and, if requested, the
#'   confidence limits. The family-wise test is stored in `maxt`, and the
#'   reparametrised models (one per reference level) in `models`.
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-based multiple
#' testing: Examples and methods for p-value adjustment*. Wiley.
#'
#' @seealso [perm_test()] for the overall test of the factor, [perm_maxt()],
#'   [boot_ci()]
#'
#' @examples
#' set.seed(1)
#' d <- sim_blocks("within", level = "trial", n_participants = 8,
#'                 n_trials = 12, seed = 1)
#' d$stimulus <- factor(sample(c("A", "B", "C"), nrow(d), replace = TRUE))
#' d$y <- d$y + 0.8 * (d$stimulus == "C")
#' m <- lm(y ~ stimulus + participant, data = d)
#'
#' # B is small here to keep the example fast; use the default for reporting
#' perm_pairwise(m, "stimulus", unit = "trial",
#'               exchange = exch_signflip("participant"), B = 99, seed = 1)
#' @export
perm_pairwise <- function(model,
                          factor,
                          term = factor,
                          unit,
                          exchange,
                          pairs = NULL,
                          adjust = c("westfall-young", "holm", "none"),
                          time = NULL,
                          B = 4999,
                          seed = NULL,
                          cluster = NULL,
                          strata = NULL,
                          ci = c("expanded", "percentile", "bca", "basic"),
                          level = 0.95,
                          data = NULL) {
  call <- match.call()
  check_model(model)
  check_column_name(factor, "factor")
  check_column_name(term, "term")
  adjust <- match.arg(adjust)
  ci <- match.arg(ci)
  check_count(B, "B", min = 1)
  data <- model_data(model, data)

  if (!factor %in% names(data)) {
    cli::cli_abort("Column {.field {factor}} not found in the data.")
  }
  values <- data[[factor]]
  if (!is.factor(values) && !is.character(values)) {
    cli::cli_abort(c(
      "{.field {factor}} must be a factor or a character column.",
      "x" = "It is of class {.cls {class(values)}}."
    ))
  }
  lv <- if (is.factor(values)) levels(droplevels(values)) else
    sort(unique(values))
  if (length(lv) < 2) {
    cli::cli_abort("{.field {factor}} must have at least two levels.")
  }

  labels <- attr(stats::terms(model), "term.labels")
  if (!term %in% labels) {
    cli::cli_abort(c(
      "{.val {term}} is not a fixed-effect term of the model.",
      "i" = "Fixed-effect terms: {.val {labels}}."
    ))
  }
  term_vars <- label_vars(model, term)
  if (!factor %in% term_vars) {
    cli::cli_abort("{.val {term}} does not contain {.field {factor}}.")
  }
  pairs <- resolve_pairs(pairs, lv)

  # Indicator columns, one per level
  dummies <- stats::setNames(
    make.names(paste0("pw.", factor, ".", lv), unique = TRUE), lv)
  clash <- intersect(dummies, names(data))
  if (length(clash) > 0) {
    cli::cli_abort("Column{?s} {.field {clash}} already exist{?s/} in the
                    data; rename {?it/them}.")
  }
  data_pw <- data
  for (l in lv) data_pw[[dummies[[l]]]] <- as.numeric(values == l)

  # One reparametrised model per reference level
  refs <- unique(pairs[, 1])
  models <- lapply(refs, function(ref) {
    f <- reparam_formula(stats::formula(model), factor,
                         dummies[setdiff(lv, ref)])
    refit <- make_refitter(model, formula = f)
    fit <- suppressWarnings(suppressMessages(
      tryCatch(refit(data_pw), error = function(e) e)))
    if (inherits(fit, "error")) {
      cli::cli_abort(c(
        "The model with reference level {.val {ref}} could not be fitted.",
        "x" = conditionMessage(fit)
      ))
    }
    fit
  })
  names(models) <- refs

  comparisons <- paste(pairs[, 2], "-", pairs[, 1])
  specs <- lapply(seq_len(nrow(pairs)), function(i) {
    m <- models[[pairs[i, 1]]]
    target <- c(setdiff(term_vars, factor), dummies[[pairs[i, 2]]])
    lab <- Filter(function(x) setequal(label_vars(m, x), target),
                  attr(stats::terms(m), "term.labels"))
    perm_spec(m, lab[[1]], data = data_pw)
  })
  names(specs) <- comparisons

  fam <- do.call(perm_maxt, c(specs, list(
    unit = unit, exchange = exchange, method = "freedman-lane", time = time,
    adjust = "step-down", B = B, seed = seed)))
  res <- fam$results
  single <- !grepl(",", res$coef, fixed = TRUE)

  out <- data.frame(
    comparison = comparisons,
    reference = pairs[, 1],
    level = pairs[, 2],
    estimate = vapply(seq_len(nrow(pairs)), function(i) {
      if (!single[i]) return(NA_real_)
      fixed_estimates(models[[pairs[i, 1]]])[[res$coef[i]]]
    }, numeric(1)),
    stat = res$stat,
    statistic = res$statistic,
    p.value = res$p.value,
    p.adjusted = switch(adjust,
                        "westfall-young" = res$p.adjusted,
                        holm = stats::p.adjust(res$p.value, "holm"),
                        none = res$p.value),
    stringsAsFactors = FALSE
  )

  if (!is.null(cluster)) {
    if (!all(single)) {
      cli::cli_abort(c(
        "Confidence intervals are available only when each pair is compared
         through a single coefficient.",
        "i" = "{.val {term}} has several coefficients per pair; set
               {.code cluster = NULL}."
      ))
    }
    out$conf.low <- NA_real_
    out$conf.high <- NA_real_
    for (ref in refs) {
      rows <- which(pairs[, 1] == ref)
      cf <- res$coef[rows]
      bt <- boot_ci(models[[ref]], cluster, resample = "case", ci = ci,
                    level = level, strata = strata, terms = cf, B = B,
                    seed = seed, data = data_pw)
      iv <- stats::confint(bt, type = ci, level = level)[cf, , drop = FALSE]
      out$conf.low[rows] <- iv[, 1]
      out$conf.high[rows] <- iv[, 2]
    }
  }

  structure(
    list(results = out, maxt = fam, models = models, factor = factor,
         term = term,
         adjust = adjust, combine = fam$combine, ci = if (!is.null(cluster)) ci,
         level = level, B = B, seed = seed, call = call),
    class = "lmmr_pairwise"
  )
}

# Variables of a model term label.
label_vars <- function(model, label) {
  factors <- attr(stats::terms(model), "factors")
  rownames(factors)[factors[, label] > 0]
}

# Pairs of levels as a two-column character matrix (reference, compared).
resolve_pairs <- function(pairs, lv, call = rlang::caller_env()) {
  if (is.null(pairs)) {
    return(t(utils::combn(lv, 2)))
  }
  if (is.list(pairs)) pairs <- do.call(rbind, pairs)
  if (!is.matrix(pairs) || ncol(pairs) != 2) {
    cli::cli_abort("{.arg pairs} must be a two-column matrix or a list of
                    pairs of levels.", call = call)
  }
  pairs <- matrix(as.character(pairs), ncol = 2)
  bad <- setdiff(as.vector(pairs), lv)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "Unknown level{?s} in {.arg pairs}: {.val {bad}}.",
      "i" = "Levels: {.val {lv}}."
    ), call = call)
  }
  if (any(pairs[, 1] == pairs[, 2])) {
    cli::cli_abort("Each pair in {.arg pairs} must contain two different
                    levels.", call = call)
  }
  pairs
}

# Replace a factor in a model formula (fixed and random parts) by the sum of
# its indicator variables.
reparam_formula <- function(f, factor, dummies) {
  repl <- str2lang(paste0("(", paste(dummies, collapse = " + "), ")"))
  out <- do.call("substitute", list(f, stats::setNames(list(repl), factor)))
  stats::as.formula(out, env = environment(f))
}

# Methods -------------------------------------------------------------------

#' @export
print.lmmr_pairwise <- function(x, digits = 3, ...) {
  fam <- x$maxt
  adj <- switch(x$adjust,
                "westfall-young" = paste0("Westfall-Young ", fam$combine,
                                          ", step-down"),
                holm = "Holm", none = "none")
  cat("\nPairwise permutation tests of `", x$term, "` between levels of `",
      x$factor, "`\n\n", sep = "")
  cat("Method:     ", method_label(fam$method), "\n", sep = "")
  cat("Exchange:   ", format(fam$exchange, unit = fam$unit), " (",
      fam$n_units, " units", if (!is.null(fam$exchange$block)) {
        paste0(" in ", fam$n_blocks, " blocks")
      }, ")\n", sep = "")
  cat("Adjustment: ", adj, "\n", sep = "")
  cat("Resamples:  ", fam$B_used, " of ", fam$B, " permutations used\n\n",
      sep = "")
  res <- x$results
  tab <- data.frame(Comparison = res$comparison, check.names = FALSE)
  if (any(!is.na(res$estimate))) {
    tab$Estimate <- formatC(res$estimate, digits = digits, format = "f")
  }
  if (!is.null(res$conf.low)) {
    pct <- paste0(format(100 * x$level), "% CI")
    tab[[pct]] <- paste0("[", formatC(res$conf.low, digits = digits,
                                      format = "f"), ", ",
                         formatC(res$conf.high, digits = digits,
                                 format = "f"), "]")
  }
  tab$Statistic <- paste(res$stat, "=", formatC(res$statistic,
                                                 digits = digits, format = "f"))
  tab$`p (unadjusted)` <- vapply(res$p.value, format_p, character(1))
  tab$`p (adjusted)` <- vapply(res$p.adjusted, format_p, character(1))
  print(tab, row.names = FALSE, right = FALSE)
  if (!is.null(res$conf.low)) {
    cat("\nIntervals: ", ci_label(x$ci),
        " bootstrap, not adjusted for multiplicity.\n", sep = "")
  }
  if (fam$B < 1000) {
    cat("\nNote: fewer than 1000 permutations; treat these results as",
        "exploratory.\n")
  }
  cat("\n")
  invisible(x)
}

#' Tidy pairwise permutation tests
#'
#' @param x An object of class `lmmr_pairwise`.
#' @param ... Not used.
#' @return A data frame with one row per pair.
#' @export
tidy.lmmr_pairwise <- function(x, ...) {
  out <- x$results
  out$adjust <- x$adjust
  out
}
