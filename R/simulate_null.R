# Simulation of new responses under the null hypothesis, from the user's own
# fitted model with the tested coefficient set to zero.

# Prepare everything needed to draw responses repeatedly; returns a function
# of no arguments that returns a new response vector.
make_null_simulator <- function(model, data, coef, null, ar1 = NULL,
                                series = NULL, time = NULL,
                                call = rlang::caller_env()) {
  response <- response_name(model, call = call)
  X <- stats::model.matrix(model)
  beta <- if (is_mixed(model)) lme4::fixef(model) else stats::coef(model)
  beta <- beta[colnames(X)]
  beta[is.na(beta)] <- 0
  beta[coef] <- 0
  eta_fixed <- as.vector(X %*% beta)

  sigma <- stats::sigma(model)

  # Random effects: b = Lambda u, with Lambda built from (modified) theta
  if (is_mixed(model)) {
    theta <- lme4::getME(model, "theta")
    if (null == "sharp") {
      for (cf in coef) theta <- zero_theta_rows(model, theta, cf)
    }
    Lt <- lme4::getME(model, "Lambdat")
    Lt@x <- theta[lme4::getME(model, "Lind")]
    Z <- lme4::getME(model, "Z")
    draw_re <- function() {
      u <- stats::rnorm(ncol(Z), 0, sigma)
      as.vector(Z %*% Matrix::crossprod(Lt, u))
    }
  } else {
    draw_re <- function() 0
  }

  draw_resid <- if (is.null(ar1)) {
    function() stats::rnorm(length(eta_fixed), 0, sigma)
  } else {
    ar1_resid_fun(data, series, time, ar1, sigma, call = call)
  }

  function() eta_fixed + draw_re() + draw_resid()
}

# Zero the rows of the relative covariance factor that belong to `coef`, so
# that the random slope of the tested term has zero variance and covariance.
zero_theta_rows <- function(model, theta, coef) {
  cnms <- lme4::getME(model, "cnms")
  pos <- 0L
  for (cn in cnms) {
    k <- length(cn)
    target <- match(coef, cn)
    for (c in seq_len(k)) {
      for (r in c:k) {
        pos <- pos + 1L
        if (!is.na(target) && r == target) theta[pos] <- 0
      }
    }
  }
  theta
}

# AR(1) residuals within each series, ordered by time, with marginal sd.
ar1_resid_fun <- function(data, series, time, phi, sigma,
                          call = rlang::caller_env()) {
  if (is.null(series) || is.null(time)) {
    cli::cli_abort(
      "{.arg series} and {.arg time} are required when {.arg ar1} is
       supplied.",
      call = call
    )
  }
  check_column_name(series, "series", call = call)
  check_column_name(time, "time", call = call)
  missing <- setdiff(c(series, time), names(data))
  if (length(missing) > 0) {
    cli::cli_abort("Column{?s} {.field {missing}} not found in the data.",
                   call = call)
  }
  ord <- order(data[[series]], data[[time]])
  s_ord <- as.character(data[[series]])[ord]
  starts <- c(TRUE, s_ord[-1] != s_ord[-length(s_ord)])
  n <- length(ord)
  innov <- sqrt(1 - phi^2)
  function() {
    z <- stats::rnorm(n)
    e <- numeric(n)
    e[1] <- z[1]
    for (i in seq_len(n)[-1]) {
      e[i] <- if (starts[i]) z[i] else phi * e[i - 1] + innov * z[i]
    }
    out <- numeric(n)
    out[ord] <- sigma * e
    out
  }
}

response_name <- function(model, call = rlang::caller_env()) {
  lhs <- stats::formula(model)[[2]]
  if (!is.name(lhs)) {
    cli::cli_abort(c(
      "The response must be a column of the data, not an expression
       ({.code {deparse(lhs)}}).",
      "i" = "Create the transformed response as a column before fitting the
             model."
    ), call = call)
  }
  as.character(lhs)
}
