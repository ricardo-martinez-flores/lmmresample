#' Declare exchangeability for permutation tests
#'
#' These functions describe which units may be exchanged when building the
#' permutation distribution. They are passed to the `exchange` argument of
#' [perm_test()], together with `unit`, the column identifying the units
#' whose labels are permuted.
#'
#' * `exch_within()`: units are exchanged only within each level of `block`.
#'   Use it for within-participant designs, where trials (`unit`) are
#'   permuted within each participant (`block`).
#' * `exch_free()`: units are exchanged freely across the whole data set, or
#'   within each level of `strata` if supplied. Use it for between-group
#'   designs, where participants (`unit`) are permuted across groups.
#' * `exch_signflip()`: the two levels of the tested variable are swapped for
#'   all units of a randomly chosen subset of blocks, which reverses the sign
#'   of each selected participant's effect. Use it in within-participant
#'   designs when participants may differ in the size of their effect (random
#'   slopes) and the hypothesis concerns the mean effect. It assumes that
#'   individual effects are symmetrically distributed around the mean under
#'   the null hypothesis.
#'
#' `exch_within()` is exact when the tested variable has no effect in any
#' participant (sharp null). When participants differ in their effect, it can
#' be liberal for the mean effect; `exch_signflip()` is then preferable. Use
#' [perm_calibrate()] to check either choice for a given design.
#'
#' Each unit must belong to a single block (or stratum), and the permuted
#' variable must be constant within each unit.
#'
#' @param block Name of the column defining the blocks within which units are
#'   exchanged (or whose labels are swapped), typically the participant
#'   identifier.
#' @param strata Optional name of a column defining strata within which units
#'   are exchanged (e.g. recording site). `NULL` (the default) exchanges units
#'   freely.
#'
#' @return An object of class `lmmr_exchange`.
#'
#' @examples
#' # Within-participant design: permute trials within participants
#' exch_within("participant")
#'
#' # Between-group design: permute participants freely
#' exch_free()
#'
#' # Between-group design stratified by site
#' exch_free(strata = "site")
#'
#' # Within-participant design, inference on the mean effect when
#' # participants differ in their effect
#' exch_signflip("participant")
#' @name exchange
NULL

#' @rdname exchange
#' @export
exch_within <- function(block) {
  check_column_name(block, "block")
  structure(list(type = "within", block = block), class = "lmmr_exchange")
}

#' @rdname exchange
#' @export
exch_signflip <- function(block) {
  check_column_name(block, "block")
  structure(list(type = "signflip", block = block), class = "lmmr_exchange")
}

#' @rdname exchange
#' @export
exch_free <- function(strata = NULL) {
  if (!is.null(strata)) check_column_name(strata, "strata")
  structure(list(type = "free", block = strata), class = "lmmr_exchange")
}

#' @export
format.lmmr_exchange <- function(x, unit = NULL, ...) {
  what <- if (is.null(unit)) "units" else paste0("`", unit, "`")
  if (x$type == "within") {
    paste0(what, " permuted within `", x$block, "`")
  } else if (x$type == "signflip") {
    paste0("labels of ", what, " swapped (sign-flipped) by `", x$block, "`")
  } else if (is.null(x$block)) {
    paste0(what, " permuted freely")
  } else {
    paste0(what, " permuted freely within strata `", x$block, "`")
  }
}

#' @export
print.lmmr_exchange <- function(x, ...) {
  cat("<exchangeability> ", format(x), "\n", sep = "")
  invisible(x)
}

# Unit table --------------------------------------------------------------

# Build one row per unit with the value of the permuted variable and its
# block, after validating the design. Returns a list with the unit table and
# a vector mapping each data row to its unit.
build_units <- function(data, term, unit, exchange, call = rlang::caller_env()) {
  if (!inherits(exchange, "lmmr_exchange")) {
    cli::cli_abort(
      "{.arg exchange} must be created with {.fn exch_within},
       {.fn exch_signflip} or {.fn exch_free}.",
      call = call
    )
  }
  check_column_name(unit, "unit", call = call)
  needed <- c(term, unit, exchange$block)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    cli::cli_abort(
      "Column{?s} {.field {missing}} not found in the data.",
      call = call
    )
  }
  for (col in needed) {
    if (anyNA(data[[col]])) {
      cli::cli_abort("Column {.field {col}} contains missing values.",
                     call = call)
    }
  }

  unit_id <- as.character(data[[unit]])
  row_unit <- match(unit_id, unique(unit_id))
  first <- !duplicated(row_unit)

  n_values <- tapply(data[[term]], row_unit, function(x) length(unique(x)))
  if (any(n_values > 1)) {
    bad <- sum(n_values > 1)
    cli::cli_abort(c(
      "{.field {term}} must be constant within each {.field {unit}}.",
      "x" = "{bad} unit{?s} contain{?s/} more than one value of
             {.field {term}}.",
      "i" = "Choose a {.arg unit} at which {.field {term}} is defined, e.g.
             the trial for a condition or the participant for a group."
    ), call = call)
  }

  units <- data.frame(
    unit = unit_id[first],
    value_row = which(first),
    stringsAsFactors = FALSE
  )
  if (!is.null(exchange$block)) {
    blk <- exchange$block
    n_blocks <- tapply(as.character(data[[blk]]), row_unit,
                       function(x) length(unique(x)))
    if (any(n_blocks > 1)) {
      label <- if (exchange$type == "free") "stratum" else "block"
      cli::cli_abort(c(
        "Each {.field {unit}} must belong to a single {label}
         ({.field {blk}}).",
        "x" = "{sum(n_blocks > 1)} unit{?s} span{?s/} several levels of
               {.field {blk}}."
      ), call = call)
    }
    units$block <- as.character(data[[blk]])[first]
  } else {
    units$block <- "all"
  }

  list(units = units, row_unit = row_unit)
}

# Number of distinct relabellings of the units (log scale), and number of
# blocks in which the variable has a single value (no information).
count_permutations <- function(values, block, type = "within") {
  groups <- split(values, block)
  uninformative <- sum(vapply(groups, function(v) {
    length(unique(as.character(v))) < 2
  }, logical(1)))
  log_n <- if (type == "signflip") {
    (length(groups) - 1) * log(2)
  } else {
    sum(vapply(groups, function(v) {
      counts <- table(as.character(v))
      lfactorial(length(v)) - sum(lfactorial(counts))
    }, numeric(1)))
  }
  list(log_n = log_n, uninformative = uninformative, n_blocks = length(groups))
}

# Matrix of relabelled values, coded as indices into `levels(codes)`: column
# b gives the value each unit receives in resample b.
generate_relabels <- function(codes, block, type, B) {
  if (type == "signflip") {
    blocks <- unique(block)
    idx <- match(block, blocks)
    out <- matrix(codes, nrow = length(codes), ncol = B)
    flips <- matrix(stats::runif(length(blocks) * B) < 0.5,
                    nrow = length(blocks))
    swap <- flips[idx, , drop = FALSE]
    out[swap] <- 3L - out[swap]
    return(out)
  }
  perms <- generate_permutations(block, B)
  matrix(codes[perms], nrow = length(codes))
}

# Matrix of permuted unit indices: column b gives, for each unit, the unit
# whose value it receives in permutation b.
generate_permutations <- function(block, B) {
  idx_by_block <- split(seq_along(block), block)
  perms <- matrix(seq_along(block), nrow = length(block), ncol = B)
  for (idx in idx_by_block) {
    if (length(idx) > 1) {
      perms[idx, ] <- vapply(seq_len(B), function(b) idx[sample.int(length(idx))],
                             integer(length(idx)))
    }
  }
  perms
}

# Encode unit values as integer codes into a vector of distinct values.
encode_values <- function(values, exchange, call = rlang::caller_env()) {
  distinct <- unique(values)
  if (is.factor(distinct)) distinct <- distinct[order(as.integer(distinct))]
  if (exchange$type == "signflip" && length(distinct) != 2) {
    cli::cli_abort(c(
      "{.fn exch_signflip} requires a variable with exactly two values.",
      "x" = "The tested variable has {length(distinct)} distinct values."
    ), call = call)
  }
  list(codes = match(values, distinct), distinct = distinct)
}

check_column_name <- function(x, arg, call = rlang::caller_env()) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    cli::cli_abort("{.arg {arg}} must be a single column name.", call = call)
  }
  invisible(x)
}
