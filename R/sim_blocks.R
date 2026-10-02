#' Simulate block-structured repeated measurements
#'
#' Generates data with a known effect structure in which densely sampled
#' series (e.g. pupil size over time) are recorded in trials nested within
#' participants. The data can be returned at three levels of aggregation: the
#' full time series, one row per trial, or one row per participant (and
#' condition). Used in examples, tests and vignettes, and useful for power
#' analysis and for checking an analysis pipeline before applying it to real
#' data.
#'
#' @section Generative model:
#' For participant \eqn{i}, trial \eqn{j} and sample \eqn{t},
#' \deqn{y_{ijt} = s(t) + \beta_c c_{ij} + \beta_g g_i + \beta_{cg} c_{ij} g_i
#'   + b_{0i} + b_{1i} c_{ij} + u_{ij} + e_{ijt},}
#' where \eqn{s(t)} is a smooth response curve common to all trials,
#' \eqn{c_{ij}} and \eqn{g_i} are 0/1 indicators for the second level of
#' `condition` and `group`, \eqn{b_{0i} \sim N(0, \sigma_p^2)} and
#' \eqn{b_{1i} \sim N(0, \sigma_s^2)} are participant random intercepts and
#' condition slopes, \eqn{u_{ij} \sim N(0, \sigma_t^2)} is a trial random
#' intercept, and \eqn{e_{ijt}} follows a stationary AR(1) process within each
#' trial with coefficient `ar1` and marginal standard deviation
#' `sd_residual`.
#'
#' Effects are constant over time, so in a model containing `condition`
#' (and/or `group`) the corresponding fixed-effect coefficient estimates
#' `effect_condition` (and/or `effect_group`) directly.
#'
#' @param design Design of the study. `"within"`: every participant completes
#'   trials in both conditions. `"between"`: participants belong to one of two
#'   groups and there is no condition factor. `"mixed"`: both.
#' @param level Level of aggregation of the returned data. `"timeseries"`: one
#'   row per sample. `"trial"`: one row per trial, with the mean (`y`) and the
#'   maximum (`peak`) of each series. `"participant"`: trial-level values
#'   averaged within each participant (and condition, if present).
#' @param n_participants Number of participants. In `"between"` and `"mixed"`
#'   designs they are split as evenly as possible between the two groups.
#' @param n_trials Number of trials per participant: a single number, or a
#'   vector with one value per participant to simulate participants with
#'   different numbers of trials (e.g. after trial rejection). In `"within"` and
#'   `"mixed"` designs, a proportion `prop_condition` of them is assigned to
#'   condition `"B"` in random order; with balanced conditions (the default)
#'   `n_trials` must be even.
#' @param n_time Number of samples per trial.
#' @param sampling_rate Sampling rate in Hz, used to express `time` in
#'   seconds.
#' @param effect_condition Difference between the second and first level of
#'   `condition` (`"B"` minus `"A"`).
#' @param effect_group Difference between the second and first level of
#'   `group` (`"treatment"` minus `"control"`).
#' @param effect_interaction Additional effect of condition `"B"` in group
#'   `"treatment"` (mixed designs only).
#' @param sd_participant Standard deviation of participant random intercepts.
#' @param sd_slope Standard deviation of participant random slopes for
#'   `condition`. Values greater than zero mean that participants differ in
#'   their condition effect.
#' @param sd_trial Standard deviation of trial random intercepts.
#' @param sd_residual Marginal standard deviation of the within-trial
#'   residuals.
#' @param ar1 Lag-1 autocorrelation of the within-trial residuals, in
#'   \eqn{[0, 1)}.
#' @param prop_condition Proportion of each participant's trials assigned to
#'   condition `"B"` (`"within"` and `"mixed"` designs). The default, 0.5,
#'   gives balanced conditions; e.g. 0.2 mimics an oddball design with 20%
#'   targets. The number of `"B"` trials is rounded to the nearest integer and
#'   kept between 1 and `n_trials - 1`.
#' @param residual_df Degrees of freedom of a Student t distribution for the
#'   innovations of the within-trial residuals. `Inf` (the default) gives
#'   normal residuals; small values (e.g. 4 or 5) give heavy-tailed
#'   residuals. Residuals are rescaled so that their marginal standard
#'   deviation is `sd_residual`.
#' @param seed Optional integer seed. If supplied, the simulation is
#'   reproducible and the global random number generator state is left
#'   unchanged.
#'
#' @return A data frame. Columns depend on `design` and `level`:
#'   * `participant`: participant identifier (factor).
#'   * `group`: `"control"` or `"treatment"` (`"between"` and `"mixed"`).
#'   * `trial`: unique trial identifier (factor; `"timeseries"` and
#'     `"trial"` levels).
#'   * `trial_index`: position of the trial within the participant's session.
#'   * `condition`: `"A"` or `"B"` (`"within"` and `"mixed"`).
#'   * `time`: time from trial onset in seconds (`"timeseries"` level).
#'   * `y`: the signal (`"timeseries"`) or its mean over the trial.
#'   * `peak`: maximum of the series within the trial (`"trial"` and
#'     `"participant"` levels).
#'
#'   The simulation parameters are stored in the attribute `"params"`.
#'
#' @examples
#' # Full time series, within-participant design
#' d <- sim_blocks("within", n_participants = 10, n_trials = 20,
#'                 effect_condition = 0.3, seed = 1)
#' head(d)
#'
#' # One row per trial
#' d_trial <- sim_blocks("within", level = "trial", seed = 1)
#'
#' # One row per participant, between-group design
#' d_part <- sim_blocks("between", level = "participant",
#'                      effect_group = 0.5, seed = 1)
#' aggregate(y ~ group, data = d_part, FUN = mean)
#' @export
sim_blocks <- function(design = c("within", "between", "mixed"),
                       level = c("timeseries", "trial", "participant"),
                       n_participants = 30,
                       n_trials = 40,
                       n_time = 100,
                       sampling_rate = 50,
                       effect_condition = 0,
                       effect_group = 0,
                       effect_interaction = 0,
                       sd_participant = 1,
                       sd_slope = 0,
                       sd_trial = 0.5,
                       sd_residual = 1,
                       ar1 = 0.9,
                       prop_condition = 0.5,
                       residual_df = Inf,
                       seed = NULL) {
  design <- match.arg(design)
  level <- match.arg(level)

  check_count(n_participants, "n_participants", min = 2)
  if (!is.numeric(n_trials) || anyNA(n_trials) ||
      any(n_trials != round(n_trials)) || any(n_trials < 1)) {
    cli::cli_abort("{.arg n_trials} must contain whole numbers of at least 1.")
  }
  if (length(n_trials) != 1 && length(n_trials) != n_participants) {
    cli::cli_abort(c(
      "{.arg n_trials} must have length 1 or {.arg n_participants}
       ({n_participants}).",
      "x" = "It has length {length(n_trials)}."
    ))
  }
  check_count(n_time, "n_time", min = 1)
  check_number(sampling_rate, "sampling_rate", min = 0, min_inclusive = FALSE)
  check_number(effect_condition, "effect_condition")
  check_number(effect_group, "effect_group")
  check_number(effect_interaction, "effect_interaction")
  check_number(sd_participant, "sd_participant", min = 0)
  check_number(sd_slope, "sd_slope", min = 0)
  check_number(sd_trial, "sd_trial", min = 0)
  check_number(sd_residual, "sd_residual", min = 0)
  check_number(ar1, "ar1", min = 0)
  if (ar1 >= 1) {
    cli::cli_abort("{.arg ar1} must be smaller than 1, not {ar1}.")
  }
  check_number(prop_condition, "prop_condition", min = 0, max = 1,
               min_inclusive = FALSE)
  if (prop_condition >= 1) {
    cli::cli_abort("{.arg prop_condition} must be smaller than 1.")
  }
  if (!is.numeric(residual_df) || length(residual_df) != 1 ||
      is.na(residual_df) || residual_df <= 2) {
    cli::cli_abort("{.arg residual_df} must be a single number greater than 2
                    (use {.code Inf} for normal residuals).")
  }

  has_condition <- design %in% c("within", "mixed")
  has_group <- design %in% c("between", "mixed")

  if (has_condition && any(n_trials < 2)) {
    cli::cli_abort(
      "{.arg n_trials} must be at least 2 in a {.val {design}} design."
    )
  }
  if (has_condition && prop_condition == 0.5 && any(n_trials %% 2 != 0)) {
    cli::cli_abort(c(
      "{.arg n_trials} must be even in a {.val {design}} design with balanced
       conditions.",
      "i" = "Half of the trials are assigned to each condition."
    ))
  }
  if (!has_condition && (effect_condition != 0 || sd_slope != 0)) {
    cli::cli_warn(
      "{.arg effect_condition} and {.arg sd_slope} are ignored in a
       {.val between} design, which has no condition factor."
    )
  }
  if (!has_group && effect_group != 0) {
    cli::cli_warn(
      "{.arg effect_group} is ignored in a {.val within} design, which has no
       group factor."
    )
  }
  if (design != "mixed" && effect_interaction != 0) {
    cli::cli_warn(
      "{.arg effect_interaction} is only used in a {.val mixed} design."
    )
  }

  params <- list(
    design = design, level = level, n_participants = n_participants,
    n_trials = n_trials, n_time = n_time, sampling_rate = sampling_rate,
    effect_condition = effect_condition, effect_group = effect_group,
    effect_interaction = effect_interaction, sd_participant = sd_participant,
    sd_slope = sd_slope, sd_trial = sd_trial, sd_residual = sd_residual,
    ar1 = ar1, prop_condition = prop_condition, residual_df = residual_df,
    seed = seed
  )

  data <- with_seed(seed, simulate_series(params, has_condition, has_group))
  data <- aggregate_level(data, level, has_condition, has_group)
  rownames(data) <- NULL
  attr(data, "params") <- params
  data
}

# Generate the full time series -------------------------------------------

simulate_series <- function(p, has_condition, has_group) {
  n_p <- p$n_participants
  k <- rep_len(p$n_trials, n_p)
  n_tr <- sum(k)
  n_t <- p$n_time

  pid_width <- nchar(n_p)
  tid_width <- nchar(max(k))
  participant <- sprintf(paste0("p%0", pid_width, "d"), seq_len(n_p))

  # Participant-level quantities
  group01 <- if (has_group) sample(rep_len(0:1, n_p)) else rep(0L, n_p)
  b0 <- stats::rnorm(n_p, 0, p$sd_participant)
  b1 <- if (has_condition) stats::rnorm(n_p, 0, p$sd_slope) else rep(0, n_p)

  # Trial-level quantities
  trial_p <- rep(seq_len(n_p), times = k)
  trial_index <- sequence(k)
  cond01 <- if (has_condition) {
    unlist(lapply(k, function(n_k) {
      n_b <- min(max(round(n_k * p$prop_condition), 1), n_k - 1)
      sample(rep(0:1, c(n_k - n_b, n_b)))
    }))
  } else {
    rep(0L, n_tr)
  }
  u <- stats::rnorm(n_tr, 0, p$sd_trial)

  trial_mean <- p$effect_condition * cond01 +
    p$effect_group * group01[trial_p] +
    p$effect_interaction * cond01 * group01[trial_p] +
    b0[trial_p] + b1[trial_p] * cond01 + u

  # Within-trial series: common response curve + AR(1) residuals
  time <- (seq_len(n_t) - 1) / p$sampling_rate
  curve <- response_curve(n_t)
  resid <- ar1_matrix(n_t, n_tr, p$ar1, p$sd_residual, p$residual_df)
  y <- as.vector(resid) + rep(curve, times = n_tr) +
    rep(trial_mean, each = n_t)

  row_trial <- rep(seq_len(n_tr), each = n_t)
  out <- data.frame(
    participant = factor(participant[trial_p[row_trial]], levels = participant),
    trial = factor(sprintf(paste0("%s_t%0", tid_width, "d"),
                           participant[trial_p], trial_index)[row_trial]),
    trial_index = trial_index[row_trial],
    time = rep(time, times = n_tr),
    y = y
  )
  if (has_group) {
    out$group <- factor(c("control", "treatment")[group01[trial_p[row_trial]] + 1],
                        levels = c("control", "treatment"))
  }
  if (has_condition) {
    out$condition <- factor(c("A", "B")[cond01[row_trial] + 1],
                            levels = c("A", "B"))
  }
  out$trial <- factor(out$trial, levels = unique(out$trial))
  order_columns(out)
}

# Smooth, pupil-like response normalised to a peak of 1 (Erlang kernel).
response_curve <- function(n_t) {
  x <- seq(0, 1, length.out = n_t)
  shape <- x^3 * exp(-6 * x)
  if (max(shape) > 0) shape / max(shape) else shape
}

# Matrix (n_t x n_series) of stationary AR(1) series with marginal sd `sd`.
# Innovations are standard normal, or Student t rescaled to unit variance.
ar1_matrix <- function(n_t, n_series, phi, sd, df = Inf) {
  z <- if (is.finite(df)) {
    stats::rt(n_t * n_series, df) * sqrt((df - 2) / df)
  } else {
    stats::rnorm(n_t * n_series)
  }
  z <- matrix(z, nrow = n_t)
  e <- z
  e[1, ] <- z[1, ]
  if (n_t > 1 && phi > 0) {
    innov <- sqrt(1 - phi^2)
    for (t in 2:n_t) e[t, ] <- phi * e[t - 1, ] + innov * z[t, ]
  }
  sd * e
}

# Aggregation -------------------------------------------------------------

aggregate_level <- function(data, level, has_condition, has_group) {
  if (level == "timeseries") return(data)

  keys <- c("participant", if (has_group) "group", "trial", "trial_index",
            if (has_condition) "condition")
  trial_id <- data$trial
  first <- !duplicated(trial_id)
  trials <- data[first, keys, drop = FALSE]
  trials$y <- as.vector(tapply(data$y, trial_id, mean))[as.integer(trials$trial)]
  trials$peak <- as.vector(tapply(data$y, trial_id, max))[as.integer(trials$trial)]
  if (level == "trial") return(trials)

  by <- c("participant", if (has_group) "group", if (has_condition) "condition")
  agg <- stats::aggregate(trials[c("y", "peak")], by = trials[by], FUN = mean)
  agg[do.call(order, unname(as.list(agg[by]))), , drop = FALSE]
}

order_columns <- function(data) {
  first <- c("participant", "group", "trial", "trial_index", "condition",
             "time", "y", "peak")
  data[c(intersect(first, names(data)), setdiff(names(data), first))]
}
