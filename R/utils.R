# Evaluate `expr` with a temporary RNG seed, restoring the global state after.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  check_count(seed, "seed", min = 0)
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  expr
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
