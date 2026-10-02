# Resampling engines --------------------------------------------------------
#
# An engine knows how to build B resampled data sets for one test. Two
# engines are available:
#
# * "relabel": the values of the tested variable are permuted (or swapped
#   by block) across whole units. Exact for a main effect without nuisance
#   terms correlated with it.
# * "freedman-lane": the reduced model (without the tested term) is fitted,
#   and its residuals are permuted across whole units, or sign-flipped by
#   block, before being added back to its fitted values (Freedman & Lane,
#   1983; Winkler et al., 2014). Valid in the presence of nuisance terms and
#   for interactions and terms with several coefficients.

# Choose the method, resolve the coefficients and the type of statistic.
resolve_test <- function(model, data, term, coef, method,
                         call = rlang::caller_env()) {
  info <- term_info(model, term, call = call)
  if (method == "auto") {
    simple <- info$is_column && !info$in_interaction && is.null(coef) &&
      length(info$coefs) == 1
    method <- if (simple || (info$is_column && !info$in_interaction &&
                             !is.null(coef))) {
      "relabel"
    } else {
      "freedman-lane"
    }
  }
  if (method == "relabel") {
    coefs <- resolve_coef(model, term, coef, call = call)
  } else {
    if (info$in_interaction) {
      cli::cli_abort(c(
        "{.val {term}} is contained in a higher-order term
         ({.val {info$higher}}).",
        "i" = "Test the highest-order interaction first; the meaning of a
               lower-order term depends on how the variables are coded when
               the interaction is in the model."
      ), call = call)
    }
    coefs <- if (is.null(coef)) info$coefs else {
      if (!all(coef %in% info$coefs)) {
        cli::cli_abort(c(
          "{.val {coef}} is not a coefficient of {.val {term}}.",
          "i" = "Available: {.val {info$coefs}}."
        ), call = call)
      }
      coef
    }
  }
  list(method = method, coefs = coefs,
       type = if (length(coefs) > 1) "chi2" else "t")
}

# Coefficients of a fixed-effect term and its position in the hierarchy.
term_info <- function(model, term, call = rlang::caller_env()) {
  tt <- stats::terms(model)
  labels <- attr(tt, "term.labels")
  factors <- attr(tt, "factors")
  if (!term %in% labels) {
    cli::cli_abort(c(
      "{.val {term}} is not a fixed-effect term of the model.",
      "i" = "Fixed-effect terms: {.val {labels}}."
    ), call = call)
  }
  vars_of <- function(lab) rownames(factors)[factors[, lab] > 0]
  my_vars <- vars_of(term)
  higher <- Filter(function(lab) {
    lab != term && all(my_vars %in% vars_of(lab))
  }, labels)

  mm <- stats::model.matrix(model)
  assign <- attr(mm, "assign")
  coefs <- colnames(mm)[assign == match(term, labels)]
  est <- if (is_mixed(model)) lme4::fixef(model) else stats::coef(model)
  coefs <- intersect(coefs, names(est)[!is.na(est)])

  list(coefs = coefs, higher = unlist(higher),
       in_interaction = length(higher) > 0,
       is_column = length(my_vars) == 1 && my_vars == term)
}

# Reduced model formula: the tested term is removed from the fixed effects
# and from any random-effects term in which it appears as a slope, because a
# random slope of the tested term would absorb the effect into the
# predicted random effects and leak it into the fitted values.
reduced_formula <- function(model, term) {
  f <- stats::formula(model)
  env <- environment(f)
  drop_label <- function(rhs, label) {
    tt <- stats::terms(stats::as.formula(paste("~", rhs)))
    labs <- setdiff(attr(tt, "term.labels"), label)
    int <- attr(tt, "intercept") == 1
    if (length(labs) == 0) return(if (int) "1" else "0")
    paste(c(if (!int) "0", labs), collapse = " + ")
  }
  lhs <- deparse(f[[2]])
  if (!is_mixed(model)) {
    rhs <- drop_label(paste(deparse(f[[3]], width.cutoff = 500L),
                            collapse = " "), term)
    return(stats::as.formula(paste(lhs, "~", rhs), env = env))
  }
  fixed <- lme4::nobars(f)
  fixed_rhs <- drop_label(paste(deparse(fixed[[3]], width.cutoff = 500L),
                                collapse = " "), term)
  bars <- vapply(lme4::findbars(f), function(bar) {
    bar_lhs <- paste(deparse(bar[[2]], width.cutoff = 500L), collapse = " ")
    group <- paste(deparse(bar[[3]], width.cutoff = 500L), collapse = " ")
    new_lhs <- drop_label(bar_lhs, term)
    # A random-effects term left without any effect is dropped entirely
    if (new_lhs == "0") return(NA_character_)
    paste0("(", new_lhs, " | ", group, ")")
  }, character(1))
  bars <- bars[!is.na(bars)]
  stats::as.formula(paste(lhs, "~", paste(c(fixed_rhs, bars),
                                          collapse = " + ")), env = env)
}

# Random-effects terms of a fitted reduced model whose grouping factor is
# constant within each exchangeability block (i.e. strictly above the
# resampled units). Only these can be kept in the fitted values used by
# Freedman-Lane: predicted random effects at or below the units would absorb
# the tested effect.
re_form_above_units <- function(fit, data, block) {
  if (!is_mixed(fit)) return(NULL)
  bars <- lme4::findbars(stats::formula(fit))
  keep <- vapply(bars, function(bar) {
    vars <- all.vars(bar[[3]])
    if (!all(vars %in% names(data))) return(FALSE)
    g <- interaction(data[vars], drop = TRUE)
    n_per_block <- tapply(as.character(g), block,
                          function(x) length(unique(x)))
    all(n_per_block == 1)
  }, logical(1))
  if (!any(keep)) return(NA)
  terms <- vapply(bars[keep], function(bar) {
    paste0("(", paste(deparse(bar, width.cutoff = 500L), collapse = " "), ")")
  }, character(1))
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

# Build an engine for one data set -----------------------------------------

build_engine <- function(model, data, term, test, unit, exchange,
                         time = NULL, call = rlang::caller_env()) {
  if (test$method == "relabel") {
    relabel_engine(data, term, unit, exchange, call = call)
  } else {
    fl_engine(model, data, term, unit, exchange, time = time, call = call)
  }
}

relabel_engine <- function(data, term, unit, exchange,
                           call = rlang::caller_env()) {
  ub <- build_units(data, term, unit, exchange, call = call)
  values <- data[[term]][ub$units$value_row]
  enc <- encode_values(values, exchange, call = call)
  space <- count_permutations(values, ub$units$block, exchange$type)
  unit_ids <- ub$units$unit
  list(
    space = space,
    n_units = nrow(ub$units),
    unit_ids = unit_ids,
    unit_blocks = ub$units$block,
    unit_values = as.character(values),
    # Express draws of another engine (`ref`) in this engine's ordering of
    # units and coding of values.
    translate = function(draws, ref) {
      map <- match(unit_ids, ref$unit_ids)
      code_map <- match(as.character(ref$distinct), as.character(enc$distinct))
      matrix(code_map[draws[map, , drop = FALSE]], nrow = length(unit_ids))
    },
    distinct = enc$distinct,
    draw = function(B) {
      generate_relabels(enc$codes, ub$units$block, exchange$type, B)
    },
    make_data = function(draws, b) {
      newdata <- data
      newdata[[term]] <- enc$distinct[draws[, b]][ub$row_unit]
      newdata
    }
  )
}

fl_engine <- function(model, data, term, unit, exchange, time = NULL,
                      call = rlang::caller_env()) {
  response <- response_name(model, call = call)
  ub <- build_units(data, NULL, unit, exchange, call = call)
  units <- ub$units
  row_unit <- ub$row_unit

  # Reduced model fitted to these data
  reduced <- make_refitter(model, formula = reduced_formula(model, term),
                           call = call)
  fit_r <- suppressWarnings(suppressMessages(
    tryCatch(reduced(data), error = function(e) e)
  ))
  if (inherits(fit_r, "error")) {
    cli::cli_abort(c(
      "The reduced model (without {.val {term}}) could not be fitted.",
      "x" = conditionMessage(fit_r)
    ), call = call)
  }
  # Fitted values keep only random effects above the resampled units
  block_of_unit_row <- units$block[row_unit]
  fitted_r <- if (is_mixed(fit_r)) {
    re_form <- re_form_above_units(fit_r, data, block_of_unit_row)
    as.vector(stats::predict(fit_r, re.form = re_form))
  } else {
    as.vector(stats::fitted(fit_r))
  }
  resid_r <- data[[response]] - fitted_r

  # Position of each row within its unit (ordered by `time` if supplied)
  sizes <- tabulate(row_unit, nbins = nrow(units))
  if (exchange$type != "signflip") {
    if (any(sizes > 1) && is.null(time)) {
      cli::cli_abort(c(
        "Permuting residuals between units with several rows requires
         {.arg time}, so that samples are matched by their time within the
         unit.",
        "i" = "Supply the column giving the order of samples within units
               (e.g. {.code time = \"time\"}), or use {.fn exch_signflip},
               which does not move residuals between units."
      ), call = call)
    }
    same_size <- tapply(sizes, units$block, function(x) length(unique(x)) == 1)
    if (!all(same_size)) {
      cli::cli_abort(c(
        "Permuting residuals between units requires units of equal size
         within each {if (exchange$type == 'free') 'stratum' else 'block'}.",
        "x" = "Units of {.field {unit}} have different numbers of rows.",
        "i" = "Use {.fn exch_signflip}, which flips whole blocks and accepts
               units of any size, or make units the same length."
      ), call = call)
    }
  }
  tvals <- if (is.null(time)) seq_along(row_unit) else {
    check_column_name(time, "time", call = call)
    if (!time %in% names(data)) {
      cli::cli_abort("Column {.field {time}} not found in the data.",
                     call = call)
    }
    data[[time]]
  }
  ord <- order(row_unit, tvals)
  offsets <- c(0, cumsum(sizes))[seq_len(nrow(units))]
  pos <- integer(length(row_unit))
  pos[ord] <- sequence(sizes)
  if (exchange$type != "signflip" && !is.null(time) && any(sizes > 1)) {
    grid <- split(tvals[ord], row_unit[ord])
    same_grid <- tapply(seq_along(grid), units$block, function(i) {
      all(vapply(grid[i], function(g) isTRUE(all.equal(g, grid[[i[1]]])),
                 logical(1)))
    })
    if (!all(same_grid)) {
      cli::cli_abort(c(
        "Units exchanged within the same block have different time points.",
        "i" = "Residuals can only be permuted between units sampled at the
               same times; use {.fn exch_signflip} otherwise."
      ), call = call)
    }
  }

  block_of_row <- match(units$block[row_unit], unique(units$block))
  space <- count_permutations(units$unit, units$block, exchange$type)
  if (exchange$type == "signflip") space$uninformative <- 0L
  unit_ids <- units$unit
  block_ids <- unique(units$block)

  list(
    space = space,
    n_units = nrow(units),
    unit_ids = unit_ids,
    unit_blocks = units$block,
    block_ids = block_ids,
    # Express draws generated by another engine (`ref`) on the same units in
    # this engine's own ordering of units and blocks.
    translate = function(draws, ref) {
      if (exchange$type == "signflip") {
        draws[match(block_ids, ref$block_ids), , drop = FALSE]
      } else {
        map <- match(unit_ids, ref$unit_ids)
        inv <- match(ref$unit_ids, unit_ids)
        matrix(inv[draws[map, , drop = FALSE]], nrow = length(unit_ids))
      }
    },
    draw = function(B) {
      if (exchange$type == "signflip") {
        n_blocks <- length(block_ids)
        matrix(sample(c(-1L, 1L), n_blocks * B, replace = TRUE),
               nrow = n_blocks)
      } else {
        generate_permutations(units$block, B)
      }
    },
    make_data = function(draws, b) {
      newdata <- data
      if (exchange$type == "signflip") {
        e <- resid_r * draws[block_of_row, b]
      } else {
        source_unit <- draws[row_unit, b]
        e <- resid_r[ord[offsets[source_unit] + pos]]
      }
      newdata[[response]] <- fitted_r + e
      newdata
    }
  )
}
