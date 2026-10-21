#' Bootstrap confidence intervals for fixed effects
#'
#' Computes bootstrap confidence intervals for the fixed-effect coefficients
#' of a linear model by resampling whole clusters (typically participants),
#' so that any dependence within clusters, such as autocorrelated samples
#' within trials, is preserved.
#'
#' @details
#' Two resampling schemes are available:
#'
#' * `resample = "case"`: clusters are drawn with replacement and the model
#'   is refitted to the resampled data. Clusters drawn more than once are
#'   given distinct identifiers, as are grouping factors nested within them
#'   (e.g. trials), so that each copy is treated as a different participant.
#' * `resample = "wild"`: the model is refitted to
#'   \eqn{y^* = X\hat\beta + v_c (y - X\hat\beta)}, where \eqn{X\hat\beta} are
#'   the fitted values of the fixed effects and \eqn{v_c = \pm 1} is drawn
#'   independently for each cluster (Rademacher weights). It keeps the
#'   design fixed, is robust to heteroscedasticity and behaves better than
#'   the case bootstrap with few clusters.
#'
#' Four interval types are available, and all replicates are stored, so
#' [confint()] can return another type without refitting (BCa only if the
#' jackknife was computed):
#'
#' * `"expanded"` (default): the percentile interval taken at a wider nominal
#'   level, so that it matches a *t* interval with \eqn{G - 1} degrees of
#'   freedom, where \eqn{G} is the number of clusters (expanded percentile
#'   interval; Hesterberg, 2015). Plain percentile intervals are too narrow
#'   with few clusters; the correction removes most of that shortfall and has
#'   little effect with many clusters (about 2% wider with 30 clusters).
#'   The limits are the quantiles at
#'   \eqn{\Phi(\mp \sqrt{G / (G - 1)}\, t_{1 - \alpha/2,\, G - 1})}.
#' * `"percentile"`: quantiles of the bootstrap distribution.
#' * `"bca"`: bias-corrected and accelerated (Efron, 1987), with acceleration
#'   estimated by a leave-one-cluster-out jackknife. It corrects for bias and
#'   skewness and is recommended with at least 2000 replicates. The jackknife
#'   requires one extra refit per cluster and is only computed when
#'   `ci = "bca"`.
#' * `"basic"`: reflection of the percentile interval around the estimate.
#'
#' With `strata`, case resampling draws clusters separately within each
#' stratum, so every bootstrap sample keeps the observed number of clusters
#' per stratum (e.g. 20 controls and 10 patients). Use it when clusters belong
#' to groups that are compared, such as patients and controls: unstratified
#' resampling would produce samples with unbalanced or empty groups. The
#' expanded percentile interval then uses \eqn{G - S} degrees of freedom,
#' where \eqn{S} is the number of strata. The wild bootstrap keeps the design
#' fixed, so strata are not needed there and `strata` is ignored.
#'
#' The case bootstrap requires clusters to be the top level of the design.
#' When the model contains random effects crossed with `cluster` (e.g. items
#' shared by all participants), the intervals reflect sampling of clusters
#' only and a warning is given; [lme4::bootMer()] provides a parametric
#' alternative.
#'
#' @param model A linear model fitted with [stats::lm()] or [lme4::lmer()].
#' @param cluster Name of the column identifying the clusters to resample.
#' @param resample Resampling scheme: `"case"` or `"wild"`. See Details.
#' @param ci Interval type: `"expanded"` (default), `"percentile"`, `"bca"`
#'   or `"basic"`. See Details.
#' @param level Confidence level.
#' @param strata Optional name of a column that is constant within each
#'   cluster (e.g. `"group"`). Case resampling is then done within its levels.
#'   See Details.
#' @param terms Names of the fixed-effect coefficients. By default all except
#'   the intercept and the coefficients of terms involving `cluster` (or
#'   grouping factors nested in it), which are not estimable in resampled
#'   data.
#' @param B Number of bootstrap replicates. The default, 4999, is suitable for
#'   reporting, in particular with BCa intervals.
#' @param seed Optional integer seed. The global random number generator
#'   state is left unchanged.
#' @param data The data used to fit `model`, if it cannot be recovered from
#'   the model call.
#'
#' @return An object of class `lmmr_boot` with methods for [print()],
#'   [confint()], [tidy()][generics::tidy] and [plot()].
#'
#' @references
#' Hesterberg, T. C. (2015). What teachers should know about the bootstrap:
#' Resampling in the undergraduate statistics curriculum. *The American
#' Statistician, 69*(4), 371-386. \doi{10.1080/00031305.2015.1089789}
#'
#' Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
#' American Statistical Association, 82*(397), 171-185.
#' \doi{10.1080/01621459.1987.10478410}
#'
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based
#' improvements for inference with clustered errors. *The Review of Economics
#' and Statistics, 90*(3), 414-427. \doi{10.1162/rest.90.3.414}
#'
#' @seealso [diag_loso()], [perm_test()]
#'
#' @examples
#' d <- sim_blocks("within", level = "trial", n_participants = 12,
#'                 n_trials = 10, effect_condition = 0.5, seed = 1)
#' m <- lm(y ~ condition + participant, data = d)
#'
#' # B is small here to keep the example fast; use the default for reporting
#' bt <- boot_ci(m, cluster = "participant", terms = "conditionB", ci = "bca",
#'               B = 99, seed = 1)
#' bt
#' confint(bt, type = "bca")
#' @export
boot_ci <- function(model,
                    cluster,
                    resample = c("case", "wild"),
                    ci = c("expanded", "percentile", "bca", "basic"),
                    level = 0.95,
                    strata = NULL,
                    terms = NULL,
                    B = 4999,
                    seed = NULL,
                    data = NULL) {
  call <- match.call()
  check_model(model)
  check_column_name(cluster, "cluster")
  resample <- match.arg(resample)
  ci <- match.arg(ci)
  check_number(level, "level", min = 0, max = 1, min_inclusive = FALSE)
  if (level >= 1) cli::cli_abort("{.arg level} must be smaller than 1.")
  check_count(B, "B", min = 2)
  data <- model_data(model, data)
  if (!cluster %in% names(data)) {
    cli::cli_abort("Column {.field {cluster}} not found in the data.")
  }
  if (anyNA(data[[cluster]])) {
    cli::cli_abort("Column {.field {cluster}} contains missing values.")
  }
  est <- fixed_estimates(model)
  n_clusters <- length(unique(data[[cluster]]))
  if (n_clusters < 3) {
    cli::cli_abort("At least 3 clusters are needed; found {n_clusters}.")
  }
  stratum_of <- cluster_strata(data, cluster, strata)
  n_strata <- if (is.null(stratum_of) || resample == "wild") 1L else
    length(unique(stratum_of))

  nested <- nested_grouping(model, data, cluster)
  if (is.null(terms)) terms <- default_terms(model, c(cluster, nested))
  terms <- check_terms(terms, est)
  refit <- make_refitter(model)
  make_data <- if (resample == "case") {
    case_resampler(data, cluster, nested, stratum_of)
  } else {
    wild_resampler(model, data, cluster)
  }

  ids <- with_seed(seed, make_data$draw(B))
  use_progress <- requireNamespace("progressr", quietly = TRUE)
  if (use_progress) p <- progressr::progressor(steps = B)
  results <- future.apply::future_lapply(seq_len(B), function(b) {
    res <- safe_refit_estimates(refit, make_data$build(ids, b), terms)
    if (use_progress) p()
    res
  }, future.seed = FALSE, future.packages = refit_packages(model))
  boot <- matrix(vapply(results, function(r) r$est, numeric(length(terms))),
                 ncol = length(terms), byrow = TRUE,
                 dimnames = list(NULL, terms))
  status <- vapply(results, function(r) r$status, character(1))
  used <- stats::complete.cases(boot)
  if (sum(!used) / B > 0.05) {
    cli::cli_warn(c(
      "{sum(!used)} of {B} bootstrap refits failed or did not converge and
       were excluded.",
      "i" = "Consider simplifying the random-effects structure."
    ))
  }

  # Jackknife (leave one cluster out) only needed for BCa acceleration
  jack <- if (ci == "bca") {
    loso_refits(refit, data, cluster, terms,
                packages = refit_packages(model))$estimates
  } else {
    NULL
  }

  out <- structure(
    list(
      estimate = est[terms],
      boot = boot,
      status = status,
      jackknife = jack,
      resample = resample,
      ci = ci,
      level = level,
      cluster = cluster,
      n_clusters = n_clusters,
      strata = if (resample == "case") strata else NULL,
      n_strata = n_strata,
      B = B,
      B_used = sum(used),
      seed = seed,
      call = call
    ),
    class = "lmmr_boot"
  )
  out$intervals <- boot_intervals(out, ci, level)
  out
}

# Grouping columns of the random effects nested within `cluster`, which must
# be relabelled when a cluster is drawn more than once. Crossed grouping
# factors trigger a warning.
nested_grouping <- function(model, data, cluster) {
  if (!is_mixed(model)) return(character(0))
  bars <- reformulas::findbars(stats::formula(model))
  cl <- as.character(data[[cluster]])
  nested <- character(0)
  for (bar in bars) {
    vars <- intersect(all.vars(bar[[3]]), names(data))
    # Groupings that include the cluster become unique once it is relabelled
    if (cluster %in% vars || length(vars) == 0) next
    g <- as.character(interaction(data[vars], drop = TRUE))
    n_per_level <- tapply(cl, g, function(x) length(unique(x)))
    if (all(n_per_level == 1)) {
      nested <- c(nested, vars)
    } else {
      label <- paste(vars, collapse = ":")
      cli::cli_warn(c(
        "Random effects grouped by {.field {label}} are crossed with
         {.field {cluster}}.",
        "i" = "Bootstrap intervals reflect sampling of {.field {cluster}}
               only, not of {.field {label}}. {.fn lme4::bootMer} provides a
               parametric alternative."
      ), call = rlang::caller_env(2))
    }
  }
  unique(nested)
}

# Stratum of each cluster (in order of first appearance), or NULL. The
# stratum column must be constant within clusters and every stratum must
# contain at least two clusters.
cluster_strata <- function(data, cluster, strata, call = rlang::caller_env()) {
  if (is.null(strata)) return(NULL)
  check_column_name(strata, "strata", call = call)
  if (!strata %in% names(data)) {
    cli::cli_abort("Column {.field {strata}} not found in the data.",
                   call = call)
  }
  if (anyNA(data[[strata]])) {
    cli::cli_abort("Column {.field {strata}} contains missing values.",
                   call = call)
  }
  cl <- as.character(data[[cluster]])
  st <- as.character(data[[strata]])
  per_cluster <- tapply(st, factor(cl, levels = unique(cl)),
                        function(x) length(unique(x)))
  if (any(per_cluster > 1)) {
    cli::cli_abort(c(
      "{.field {strata}} must be constant within each {.field {cluster}}.",
      "i" = "{sum(per_cluster > 1)} cluster{?s} belong{?s/} to more than one
             stratum."
    ), call = call)
  }
  stratum_of <- st[match(unique(cl), cl)]
  size <- table(stratum_of)
  if (any(size < 2)) {
    cli::cli_abort(c(
      "Every stratum needs at least 2 clusters.",
      "x" = "Too few in: {.val {names(size)[size < 2]}}."
    ), call = call)
  }
  stratum_of
}

case_resampler <- function(data, cluster, nested, stratum_of = NULL) {
  cl <- as.character(data[[cluster]])
  ids <- unique(cl)
  rows_of <- split(seq_len(nrow(data)), factor(cl, levels = ids))
  # Positions of the clusters of each stratum; each position is filled with a
  # cluster drawn from the same stratum
  groups <- if (is.null(stratum_of)) list(seq_along(ids)) else
    split(seq_along(ids), factor(stratum_of, levels = unique(stratum_of)))
  list(
    draw = function(B) {
      out <- matrix(0L, nrow = length(ids), ncol = B)
      for (g in groups) {
        out[g, ] <- g[sample.int(length(g), length(g) * B, replace = TRUE)]
      }
      out
    },
    build = function(draws, b) {
      pick <- draws[, b]
      rows <- unlist(rows_of[pick], use.names = FALSE)
      copy <- rep(seq_along(pick), lengths(rows_of[pick]))
      newdata <- data[rows, , drop = FALSE]
      for (v in c(cluster, nested)) {
        newdata[[v]] <- factor(paste0(as.character(newdata[[v]]), "#", copy))
      }
      rownames(newdata) <- NULL
      newdata
    }
  )
}

wild_resampler <- function(model, data, cluster) {
  response <- response_name(model)
  xb <- if (is_mixed(model)) {
    as.vector(stats::predict(model, re.form = NA))
  } else {
    as.vector(stats::fitted(model))
  }
  e <- data[[response]] - xb
  cl <- match(as.character(data[[cluster]]), unique(as.character(data[[cluster]])))
  n_cl <- max(cl)
  list(
    draw = function(B) {
      matrix(sample(c(-1L, 1L), n_cl * B, replace = TRUE), nrow = n_cl)
    },
    build = function(draws, b) {
      newdata <- data
      newdata[[response]] <- xb + draws[cl, b] * e
      newdata
    }
  )
}

# Interval computation ---------------------------------------------------------

boot_intervals <- function(x, type, level) {
  alpha <- (1 - level) / 2
  rows <- lapply(colnames(x$boot), function(cf) {
    b <- x$boot[, cf]
    b <- b[!is.na(b)]
    est <- x$estimate[[cf]]
    lim <- switch(type,
      expanded = stats::quantile(b, expanded_probs(
        x$n_clusters, alpha, if (is.null(x$n_strata)) 1L else x$n_strata),
                                 names = FALSE, type = 6),
      percentile = stats::quantile(b, c(alpha, 1 - alpha), names = FALSE,
                                   type = 6),
      basic = 2 * est - stats::quantile(b, c(1 - alpha, alpha), names = FALSE,
                                        type = 6),
      bca = bca_limits(b, est, if (is.null(x$jackknife)) NULL else
        x$jackknife[, cf], alpha)
    )
    data.frame(term = cf, estimate = est, conf.low = lim[1],
               conf.high = lim[2], type = type, level = level)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Tail probabilities of the expanded percentile interval (Hesterberg, 2015):
# the normal quantile is replaced by a t quantile with G - 1 df, inflated by
# sqrt(G / (G - 1)).
expanded_probs <- function(n_clusters, alpha, n_strata = 1L) {
  df <- n_clusters - n_strata
  if (df < 1) return(c(alpha, 1 - alpha))
  z <- sqrt(n_clusters / df) * stats::qt(1 - alpha, df)
  stats::pnorm(c(-z, z))
}

bca_limits <- function(boot, est, jack, alpha) {
  if (is.null(jack)) {
    cli::cli_abort(c(
      "BCa intervals need the leave-one-cluster-out jackknife.",
      "i" = "Run {.fn boot_ci} with {.code ci = \"bca\"}."
    ), call = rlang::caller_env(3))
  }
  jack <- jack[!is.na(jack)]
  p_below <- mean(boot < est) + 0.5 * mean(boot == est)
  if (p_below <= 0 || p_below >= 1) {
    cli::cli_warn("BCa bias correction is undefined (all bootstrap estimates
                   on one side of the estimate); returning percentile limits.")
    return(stats::quantile(boot, c(alpha, 1 - alpha), names = FALSE,
                           type = 6))
  }
  z0 <- stats::qnorm(p_below)
  d <- mean(jack) - jack
  acc <- sum(d^3) / (6 * sum(d^2)^1.5)
  if (!is.finite(acc)) acc <- 0
  z <- stats::qnorm(c(alpha, 1 - alpha))
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - acc * (z0 + z)))
  stats::quantile(boot, adj, names = FALSE, type = 6)
}

# Methods ---------------------------------------------------------------------

#' Confidence intervals from a bootstrap
#'
#' @param object An object of class `lmmr_boot`.
#' @param parm Coefficients to return; all by default.
#' @param level Confidence level.
#' @param type Interval type: `"expanded"`, `"percentile"`, `"bca"` or
#'   `"basic"`. Defaults to the type used in [boot_ci()].
#' @param ... Not used.
#' @return A matrix with one row per coefficient and the lower and upper
#'   limits, as returned by [stats::confint()].
#' @export
confint.lmmr_boot <- function(object, parm = NULL, level = object$level,
                              type = object$ci, ...) {
  type <- match.arg(type, c("expanded", "percentile", "bca", "basic"))
  iv <- boot_intervals(object, type, level)
  if (!is.null(parm)) iv <- iv[iv$term %in% parm, , drop = FALSE]
  pct <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                       trim = TRUE), " %")
  matrix(c(iv$conf.low, iv$conf.high), ncol = 2,
         dimnames = list(iv$term, pct))
}

#' @export
print.lmmr_boot <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  cat("\nBootstrap confidence intervals (", x$resample, " bootstrap of `",
      x$cluster, "`, ", x$n_clusters, " clusters",
      if (!is.null(x$strata)) paste0(" within ", x$n_strata, " strata of `",
                                     x$strata, "`"),
      ")\n\n", sep = "")
  iv <- x$intervals
  tab <- data.frame(
    Coefficient = iv$term,
    Estimate = f(iv$estimate),
    Lower = f(iv$conf.low),
    Upper = f(iv$conf.high),
    check.names = FALSE
  )
  cat(100 * x$level, "% ", ci_label(x$ci), " intervals\n", sep = "")
  print(tab, row.names = FALSE, right = FALSE)
  cat("\nReplicates: ", x$B_used, " of ", x$B, " used\n", sep = "")
  if (x$B < 1000) {
    cat("Note: fewer than 1000 replicates; treat these intervals as",
        "exploratory.\n")
  }
  cat("\n")
  invisible(x)
}

ci_label <- function(type) {
  switch(type, expanded = "expanded percentile", percentile = "percentile",
         bca = "BCa", basic = "basic")
}

#' Tidy a bootstrap
#'
#' @param x An object of class `lmmr_boot`.
#' @param type `"summary"` returns the estimates and intervals of the type
#'   used in [boot_ci()]; `"replicates"` returns all bootstrap estimates in
#'   long format.
#' @param ... Not used.
#' @return A data frame.
#' @export
tidy.lmmr_boot <- function(x, type = c("summary", "replicates"), ...) {
  type <- match.arg(type)
  if (type == "replicates") {
    return(data.frame(
      replicate = rep(seq_len(nrow(x$boot)), ncol(x$boot)),
      term = rep(colnames(x$boot), each = nrow(x$boot)),
      estimate = as.vector(x$boot)
    ))
  }
  x$intervals
}

#' Plot a bootstrap
#'
#' Shows the bootstrap distribution of each coefficient with the estimate,
#' the expanded percentile and percentile intervals, and the BCa interval if
#' it was computed, so that their differences (due to few clusters, bias or
#' skewness) can be seen.
#'
#' @param x An object of class `lmmr_boot`.
#' @param ... Not used.
#' @return A [ggplot2::ggplot()] object with one panel per coefficient.
#' @export
plot.lmmr_boot <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "to plot results.")
  reps <- tidy.lmmr_boot(x, type = "replicates")
  reps <- reps[!is.na(reps$estimate), ]
  shown <- unique(c("expanded", "percentile",
                    if (!is.null(x$jackknife)) "bca"))
  iv <- do.call(rbind, lapply(shown, function(t) boot_intervals(x, t, x$level)))
  labs <- c(expanded = "Expanded percentile", percentile = "Percentile",
            bca = "BCa")
  iv$type <- factor(labs[iv$type], levels = labs[shown])
  peak <- vapply(split(reps$estimate, reps$term), function(v) {
    if (length(unique(v)) < 2) return(1)
    max(stats::density(v)$y)
  }, numeric(1))
  iv$y <- -peak[iv$term] * 0.06 * as.integer(iv$type)
  est <- data.frame(term = names(x$estimate), estimate = unname(x$estimate))
  ggplot2::ggplot(reps, ggplot2::aes(x = .data$estimate)) +
    ggplot2::geom_density(fill = "grey85", colour = "grey40") +
    ggplot2::geom_vline(data = est, ggplot2::aes(xintercept = .data$estimate),
                        colour = "grey15") +
    ggplot2::geom_segment(
      data = iv,
      ggplot2::aes(x = .data$conf.low, xend = .data$conf.high,
                   y = .data$y, yend = .data$y, colour = .data$type),
      inherit.aes = FALSE, linewidth = 1.2
    ) +
    ggplot2::scale_colour_manual(values = c(`Expanded percentile` = "grey15",
                                            Percentile = "grey60",
                                            BCa = "#C0392B"), name = NULL) +
    ggplot2::facet_wrap(~term, scales = "free") +
    ggplot2::labs(
      x = "Bootstrap estimate", y = "Density",
      title = "Bootstrap distribution",
      subtitle = paste0(x$B_used, " ", x$resample, " bootstrap replicates; ",
                        100 * x$level, "% intervals below the axis")
    ) +
    theme_lmmr() +
    ggplot2::theme(legend.position = "bottom")
}
