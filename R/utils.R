# Evaluate `expr` with a temporary RNG seed, restoring the global state after.
with_seed <- function(seed, expr, kind = NULL) {
  if (is.null(seed)) return(expr)
  check_count(seed, "seed", min = 0)
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  old_kind <- RNGkind()
  if (!is.null(kind)) do.call(RNGkind, as.list(kind))
  on.exit({
    if (!is.null(kind)) do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  expr
}

# Series must be identified uniquely: (series, time) pairs cannot repeat.
check_series_ids <- function(data, series, time, call = rlang::caller_env()) {
  key <- paste(as.character(data[[series]]), data[[time]], sep = "\r")
  if (anyDuplicated(key)) {
    cli::cli_abort(c(
      "Some combinations of {.field {series}} and {.field {time}} occur more
       than once.",
      "i" = "Series identifiers must be unique across the data set. If trials
             are numbered within each participant, create unique identifiers,
             e.g. {.code interaction(participant, trial)}."
    ), call = call)
  }
  invisible(TRUE)
}

# Argument checks ---------------------------------------------------------

check_number <- function(x, arg, min = -Inf, max = Inf, min_inclusive = TRUE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    cli::cli_abort("{.arg {arg}} must be a single finite number.",
                   call = rlang::caller_env())
  }
  below <- if (min_inclusive) x < min else x <= min
  if (below || x > max) {
    bound <- if (min_inclusive) "at least" else "greater than"
    cli::cli_abort("{.arg {arg}} must be {bound} {min}, not {x}.",
                   call = rlang::caller_env())
  }
  invisible(x)
}

check_count <- function(x, arg, min = 0) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x)) {
    cli::cli_abort("{.arg {arg}} must be a single whole number.",
                   call = rlang::caller_env())
  }
  if (x < min) {
    cli::cli_abort("{.arg {arg}} must be at least {min}, not {x}.",
                   call = rlang::caller_env())
  }
  invisible(x)
}
