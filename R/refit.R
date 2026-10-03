# Model support ------------------------------------------------------------

check_model <- function(model, call = rlang::caller_env()) {
  ok <- (inherits(model, "lm") && !inherits(model, "glm")) ||
    inherits(model, "lmerMod")
  if (!ok) {
    cli::cli_abort(c(
      "{.arg model} must be a linear model fitted with {.fn lm} or
       {.fn lme4::lmer}.",
      "x" = "It has class {.cls {class(model)}}.",
      "i" = "Generalized linear (mixed) models are not supported."
    ), call = call)
  }
  cl <- if (inherits(model, "lmerMod")) model@call else model$call
  used <- intersect(c("weights", "offset"), names(cl))
  if (length(used) > 0 || !is.null(attr(stats::terms(model), "offset"))) {
    cli::cli_abort(c(
      "Models fitted with weights or offsets are not supported.",
      "i" = "Resampling would not carry the weights or offsets along with the
             data."
    ), call = call)
  }
  dropped <- if (inherits(model, "lmerMod")) {
    names(attr(lme4::getME(model, "X"), "col.dropped"))
  } else {
    names(stats::coef(model))[is.na(stats::coef(model))]
  }
  if (length(dropped) > 0) {
    cli::cli_abort(c(
      "The fixed-effects design of {.arg model} is rank deficient.",
      "x" = "{cli::qty(length(dropped))}Coefficient{?s} {.val {dropped}}
             could not be estimated and {?was/were} dropped.",
      "i" = "Remove redundant terms (e.g. an interaction with a variable whose
             main effect is already captured by another term) so that the
             tested coefficients are well defined."
    ), call = call)
  }
  invisible(model)
}

is_mixed <- function(model) inherits(model, "merMod")

model_call <- function(model) {
  if (is_mixed(model)) model@call else model$call
}

model_env <- function(model) {
  env <- environment(stats::formula(model))
  if (is.null(env)) globalenv() else env
}

# Recover the data the model was fitted to.
model_data <- function(model, data = NULL, call = rlang::caller_env()) {
  if (is.null(data)) {
    cl <- model_call(model)
    if (is.null(cl$data)) {
      cli::cli_abort(
        "Cannot recover the data from {.arg model}; supply it with
         {.arg data}.",
        call = call
      )
    }
    data <- tryCatch(eval(cl$data, model_env(model)), error = function(e) NULL)
    if (is.null(data)) {
      cli::cli_abort(c(
        "Cannot recover the data from {.arg model}.",
        "i" = "Supply the data used to fit the model with {.arg data}."
      ), call = call)
    }
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.", call = call)
  }
  n_model <- stats::nobs(model)
  if (nrow(data) != n_model) {
    cli::cli_abort(c(
      "The data has {nrow(data)} rows but the model was fitted to
       {n_model} observations.",
      "i" = "Remove rows with missing values before fitting the model, so
             that every row of the data is used."
    ), call = call)
  }
  data <- as.data.frame(data)
  check_alignment(model, data, call = call)
  data
}

# The data must be the rows used to fit the model, in the same order:
# several functions combine model vectors (fitted values, residuals, model
# matrix) with columns of the data.
check_alignment <- function(model, data, call = rlang::caller_env()) {
  y_model <- as.numeric(stats::model.response(stats::model.frame(model)))
  y_data <- tryCatch(
    as.numeric(eval(stats::formula(model)[[2]], data, model_env(model))),
    error = function(e) NULL
  )
  if (is.null(y_data) || length(y_data) != length(y_model) ||
      !isTRUE(all.equal(unname(y_model), unname(y_data)))) {
    cli::cli_abort(c(
      "{.arg data} does not match the data used to fit the model.",
      "i" = "Supply the same rows, in the same order, as used to fit the
             model."
    ), call = call)
  }
  invisible(TRUE)
}

# Build a function that refits the model to new data. All arguments of the
# original call except `data` are evaluated once, here, so that refits do not
# depend on the calling environment (needed for parallel workers). A
# different `formula` can be supplied (e.g. a reduced model).
make_refitter <- function(model, formula = NULL, call = rlang::caller_env()) {
  cl <- model_call(model)
  if (!is.null(cl$subset)) {
    cli::cli_abort(c(
      "Models fitted with {.arg subset} are not supported.",
      "i" = "Subset the data before fitting the model."
    ), call = call)
  }
  env <- model_env(model)
  fun <- eval(cl[[1]], env)
  args <- as.list(cl)[-1]
  args$data <- NULL
  args <- lapply(args, function(a) eval(a, env))
  args$formula <- if (is.null(formula)) stats::formula(model) else formula
  if (is_mixed(model) && isTRUE(getOption("lmmresample.fast", TRUE))) {
    # Start from the variance parameters of the original fit (same model
    # only) and skip lme4's derivative-based convergence check, which is
    # slow on large data and often raises false alarms.
    if (is.null(formula) && is.null(args$start)) {
      args$start <- list(theta = lme4::getME(model, "theta"))
    }
    if (is.null(args$control)) {
      args$control <- lme4::lmerControl(calc.derivs = FALSE)
    }
  }
  args$data <- quote(.lmmr_data)
  refit_call <- as.call(c(list(fun), args))
  function(newdata) {
    eval(refit_call, list(.lmmr_data = newdata), baseenv())
  }
}

# Packages that provide functions used in model formulas (e.g. splines for
# ns()), which must be attached in parallel workers so that refits can
# evaluate the formula.
refit_packages <- function(...) {
  models <- list(...)
  pkgs <- unlist(lapply(models, function(model) {
    f <- stats::formula(model)
    env <- model_env(model)
    fns <- unique(all.names(f))
    vapply(fns, function(fn) {
      obj <- tryCatch(get(fn, envir = env, mode = "function"),
                      error = function(e) NULL)
      if (is.null(obj) || is.primitive(obj)) return(NA_character_)
      ns <- environment(obj)
      if (is.null(ns) || !isNamespace(ns)) return(NA_character_)
      getNamespaceName(ns)
    }, character(1))
  }))
  pkgs <- unique(stats::na.omit(unname(pkgs)))
  setdiff(pkgs, c("base", "stats", "methods", "utils", "graphics",
                  "grDevices", "lmmresample"))
}

# Refit, recording convergence problems instead of failing.
# Returns list(fit = model or NULL, status = character(1)).
safe_fit <- function(refit, newdata) {
  status <- "ok"
  fit <- withCallingHandlers(
    tryCatch(refit(newdata), error = function(e) NULL),
    warning = function(w) {
      if (grepl("converge|Hessian|eigenvalue", conditionMessage(w))) {
        status <<- "nonconverged"
      }
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (grepl("singular", conditionMessage(m))) {
        if (status == "ok") status <<- "singular"
      }
      invokeRestart("muffleMessage")
    }
  )
  if (is.null(fit)) return(list(fit = NULL, status = "failed"))
  if (status == "ok" && is_mixed(fit) && lme4::isSingular(fit)) {
    status <- "singular"
  }
  list(fit = fit, status = status)
}

# Refit and extract statistics, recording problems instead of failing.
# Returns list(stat = named numeric, status = character(1)). With
# `type = "chi2"` the statistic is a single Wald chi-square for all `coefs`.
safe_refit <- function(refit, newdata, coefs, type = "t") {
  res <- safe_fit(refit, newdata)
  if (is.null(res$fit)) {
    n_stat <- if (type == "chi2") 1 else length(coefs)
    nm <- if (type == "chi2") "chi2" else coefs
    return(list(stat = stats::setNames(rep(NA_real_, n_stat), nm),
                status = "failed"))
  }
  stat <- if (type == "chi2") {
    wald_chi2(res$fit, coefs)
  } else {
    coef_stats(res$fit, coefs)
  }
  status <- res$status
  if (anyNA(stat) && status == "ok") status <- "failed"
  list(stat = stat, status = status)
}

# Wald chi-square statistic for a set of coefficients.
wald_chi2 <- function(fit, coefs) {
  est <- if (is_mixed(fit)) lme4::fixef(fit) else stats::coef(fit)
  v <- as.matrix(stats::vcov(fit))
  b <- est[coefs]
  if (anyNA(b) || !all(coefs %in% colnames(v))) return(c(chi2 = NA_real_))
  out <- tryCatch(
    as.numeric(crossprod(b, solve(v[coefs, coefs, drop = FALSE], b))),
    error = function(e) NA_real_
  )
  c(chi2 = out)
}

# Wald statistic (estimate / standard error) for the requested coefficients.
coef_stats <- function(fit, coefs) {
  est <- if (is_mixed(fit)) lme4::fixef(fit) else stats::coef(fit)
  se <- sqrt(diag(as.matrix(stats::vcov(fit))))
  names(se) <- colnames(as.matrix(stats::vcov(fit)))
  out <- est[coefs] / se[coefs]
  names(out) <- coefs
  out
}

stat_label <- function(model) "t"

# Map a variable name to its coefficient(s), refusing interactions.
resolve_coef <- function(model, term, coef = NULL, call = rlang::caller_env()) {
  tt <- stats::terms(model)
  labels <- attr(tt, "term.labels")
  factors <- attr(tt, "factors")

  if (!term %in% labels) {
    cli::cli_abort(c(
      "{.val {term}} is not a fixed-effect term of the model.",
      "i" = "Fixed-effect terms: {.val {labels}}."
    ), call = call)
  }
  if (term %in% rownames(factors)) {
    in_terms <- colnames(factors)[factors[term, ] > 0]
    interactions <- setdiff(in_terms, term)
    if (length(interactions) > 0) {
      cli::cli_abort(c(
        "{.val {term}} is involved in an interaction
         ({.val {interactions}}).",
        "i" = "Tests of terms involved in interactions require permutation of
               residuals (Freedman-Lane), planned for a future version.",
        "i" = "Fit a model without the interaction to test the main effect."
      ), call = call)
    }
  }

  mm <- stats::model.matrix(model)
  assign <- attr(mm, "assign")
  candidates <- colnames(mm)[assign == match(term, labels)]
  est <- if (is_mixed(model)) lme4::fixef(model) else stats::coef(model)
  candidates <- intersect(candidates, names(est)[!is.na(est)])

  if (is.null(coef)) {
    if (length(candidates) != 1) {
      cli::cli_abort(c(
        "{.val {term}} corresponds to {length(candidates)} coefficients
         ({.val {candidates}}).",
        "i" = "Choose one with {.arg coef}. Tests of all levels of a factor
               (F or likelihood-ratio statistics) are planned for a future
               version."
      ), call = call)
    }
    return(candidates)
  }
  if (!coef %in% candidates) {
    cli::cli_abort(c(
      "{.val {coef}} is not a coefficient of {.val {term}}.",
      "i" = "Available: {.val {candidates}}."
    ), call = call)
  }
  coef
}

# Run the refits over a matrix of permuted unit indices, in parallel through
# the future framework. `make_data(b)` returns the data for replicate b.
run_refits <- function(refit, make_data, B, coefs, type = "t",
                       packages = character(0)) {
  use_progress <- requireNamespace("progressr", quietly = TRUE)
  run <- function() {
    if (use_progress) p <- progressr::progressor(steps = B)
    future.apply::future_lapply(seq_len(B), function(b) {
      res <- safe_refit(refit, make_data(b), coefs, type)
      if (use_progress) p()
      res
    }, future.seed = FALSE, future.packages = packages)
  }
  results <- run()
  nm <- if (type == "chi2") "chi2" else coefs
  stats <- matrix(vapply(results, function(r) r$stat, numeric(length(nm))),
                  ncol = length(nm), byrow = TRUE,
                  dimnames = list(NULL, nm))
  status <- vapply(results, function(r) r$status, character(1))
  list(stats = stats, status = status)
}
