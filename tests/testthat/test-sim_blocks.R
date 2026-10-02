small <- function(...) {
  sim_blocks(n_participants = 6, n_trials = 4, n_time = 10, seed = 1, ...)
}

test_that("dimensions match the requested design and level", {
  expect_equal(nrow(small("within", "timeseries")), 6 * 4 * 10)
  expect_equal(nrow(small("within", "trial")), 6 * 4)
  expect_equal(nrow(small("within", "participant")), 6 * 2)
  expect_equal(nrow(small("between", "participant")), 6)
  expect_equal(nrow(small("mixed", "participant")), 6 * 2)
})

test_that("columns reflect the design", {
  within <- small("within")
  between <- small("between")
  mixed <- small("mixed")

  expect_true("condition" %in% names(within))
  expect_false("group" %in% names(within))
  expect_true("group" %in% names(between))
  expect_false("condition" %in% names(between))
  expect_true(all(c("group", "condition") %in% names(mixed)))

  expect_true("time" %in% names(within))
  expect_false("time" %in% names(small("within", "trial")))
  expect_true("peak" %in% names(small("within", "trial")))
})

test_that("trial identifiers are unique and condition is constant within trials", {
  d <- small("mixed", "timeseries")
  per_trial <- tapply(as.character(d$condition), d$trial,
                      function(x) length(unique(x)))
  expect_true(all(per_trial == 1))

  d_trial <- small("mixed", "trial")
  expect_false(anyDuplicated(d_trial$trial) > 0)
  expect_equal(nlevels(d_trial$trial), 6 * 4)
})

test_that("each participant has balanced conditions in within designs", {
  d <- sim_blocks("within", "trial", n_participants = 8, n_trials = 10,
                  seed = 3)
  counts <- table(d$participant, d$condition)
  expect_true(all(counts == 5))
})

test_that("groups are balanced and constant within participants", {
  d <- sim_blocks("between", "trial", n_participants = 9, n_trials = 2,
                  seed = 3)
  per_participant <- tapply(as.character(d$group), d$participant,
                            function(x) length(unique(x)))
  expect_true(all(per_participant == 1))
  sizes <- table(unique(d[c("participant", "group")])$group)
  expect_true(all(sizes %in% c(4, 5)))
})

test_that("a seed makes the simulation reproducible", {
  expect_identical(small("mixed"), small("mixed"))
  d1 <- sim_blocks("within", "trial", seed = 1)
  d2 <- sim_blocks("within", "trial", seed = 2)
  expect_false(identical(d1$y, d2$y))
})

test_that("a seed does not change the global random number generator", {
  set.seed(42)
  expected <- stats::runif(1)
  set.seed(42)
  small("within")
  expect_identical(stats::runif(1), expected)
})

test_that("the condition effect is recovered", {
  d <- sim_blocks("within", "participant", n_participants = 400,
                  n_trials = 20, n_time = 20, effect_condition = 0.5,
                  seed = 10)
  diffs <- d$y[d$condition == "B"] - d$y[d$condition == "A"]
  expect_equal(mean(diffs), 0.5, tolerance = 0.1)
})

test_that("the group effect is recovered", {
  d <- sim_blocks("between", "participant", n_participants = 2000,
                  n_trials = 4, n_time = 10, effect_group = -0.8, seed = 10)
  means <- tapply(d$y, d$group, mean)
  expect_equal(unname(means["treatment"] - means["control"]), -0.8,
               tolerance = 0.1)
})

test_that("the interaction effect is recovered", {
  d <- sim_blocks("mixed", "participant", n_participants = 2000,
                  n_trials = 4, n_time = 10, effect_interaction = 1,
                  sd_trial = 0, seed = 11)
  wide <- stats::reshape(d[c("participant", "group", "condition", "y")],
                         idvar = c("participant", "group"),
                         timevar = "condition", direction = "wide")
  diff <- wide$y.B - wide$y.A
  contrast <- mean(diff[wide$group == "treatment"]) -
    mean(diff[wide$group == "control"])
  expect_equal(contrast, 1, tolerance = 0.1)
})

test_that("residuals follow the requested AR(1) structure", {
  d <- sim_blocks("within", "timeseries", n_participants = 4, n_trials = 50,
                  n_time = 200, sd_participant = 0, sd_trial = 0, ar1 = 0.7,
                  seed = 5)
  curve <- tapply(d$y, d$time, mean)
  resid <- d$y - curve[as.character(d$time)]
  lag1 <- vapply(split(resid, d$trial), function(r) {
    stats::cor(r[-1], r[-length(r)])
  }, numeric(1))
  expect_equal(mean(lag1), 0.7, tolerance = 0.05)
})

test_that("random slopes produce between-participant variation in effects", {
  d <- sim_blocks("within", "participant", n_participants = 300,
                  n_trials = 40, n_time = 10, sd_slope = 1, sd_trial = 0,
                  sd_residual = 0.1, seed = 6)
  diffs <- d$y[d$condition == "B"] - d$y[d$condition == "A"]
  expect_equal(stats::sd(diffs), 1, tolerance = 0.15)
})

test_that("parameters are stored as an attribute", {
  d <- small("within", effect_condition = 0.2)
  params <- attr(d, "params")
  expect_equal(params$effect_condition, 0.2)
  expect_equal(params$design, "within")
})

test_that("invalid arguments give informative errors", {
  expect_error(sim_blocks("within", n_trials = 5), "must be even")
  expect_error(sim_blocks(n_participants = 1), "at least 2")
  expect_error(sim_blocks(n_participants = 2.5), "whole number")
  expect_error(sim_blocks(ar1 = 1), "smaller than 1")
  expect_error(sim_blocks(sd_residual = -1), "at least 0")
  expect_error(sim_blocks(sampling_rate = 0), "greater than 0")
  expect_error(sim_blocks(design = "crossed"))
})

test_that("ignored effects trigger warnings", {
  expect_warning(small("between", effect_condition = 1), "ignored")
  expect_warning(small("within", effect_group = 1), "ignored")
  expect_warning(small("within", effect_interaction = 1), "mixed")
})

test_that("conditions can be unbalanced", {
  d <- sim_blocks("within", "trial", n_participants = 5, n_trials = 50,
                  prop_condition = 0.2, seed = 1)
  counts <- table(d$participant, d$condition)
  expect_true(all(counts[, "B"] == 10))
  expect_true(all(counts[, "A"] == 40))

  d_odd <- sim_blocks("within", "trial", n_participants = 3, n_trials = 7,
                      prop_condition = 0.3, seed = 1)
  expect_true(all(table(d_odd$participant, d_odd$condition)[, "B"] == 2))

  d_min <- sim_blocks("within", "trial", n_participants = 3, n_trials = 10,
                      prop_condition = 0.01, seed = 1)
  expect_true(all(table(d_min$participant, d_min$condition)[, "B"] == 1))

  expect_error(sim_blocks(prop_condition = 1), "smaller than 1")
  expect_error(sim_blocks(prop_condition = 0), "greater than 0")
})

test_that("heavy-tailed residuals keep the requested standard deviation", {
  d <- sim_blocks("within", "timeseries", n_participants = 4, n_trials = 100,
                  n_time = 100, sd_participant = 0, sd_trial = 0, ar1 = 0,
                  residual_df = 5, seed = 7)
  curve <- tapply(d$y, d$time, mean)
  r <- d$y - curve[as.character(d$time)]
  z <- (r - mean(r)) / stats::sd(r)
  expect_equal(stats::sd(r), 1, tolerance = 0.05)
  expect_gt(mean(z^4) - 3, 2)

  d_norm <- sim_blocks("within", "timeseries", n_participants = 4,
                       n_trials = 100, n_time = 100, sd_participant = 0,
                       sd_trial = 0, ar1 = 0, seed = 7)
  r_norm <- d_norm$y - tapply(d_norm$y, d_norm$time, mean)[
    as.character(d_norm$time)]
  z_norm <- (r_norm - mean(r_norm)) / stats::sd(r_norm)
  expect_lt(abs(mean(z_norm^4) - 3), 0.3)

  expect_error(sim_blocks(residual_df = 2), "greater than 2")
})

test_that("participants can have different numbers of trials", {
  d <- sim_blocks("within", "trial", n_participants = 3,
                  n_trials = c(6, 20, 50), prop_condition = 0.2, seed = 1)
  expect_equal(as.vector(table(d$participant)), c(6, 20, 50))
  counts_b <- as.vector(table(d$participant, d$condition)[, "B"])
  expect_equal(counts_b, c(1, 4, 10))
  expect_false(anyDuplicated(d$trial) > 0)

  ts <- sim_blocks("within", n_participants = 2, n_trials = c(4, 8),
                   n_time = 5, seed = 1)
  expect_equal(nrow(ts), (4 + 8) * 5)

  expect_error(sim_blocks(n_participants = 3, n_trials = c(10, 20)),
               "length 1 or")
  expect_error(sim_blocks(n_participants = 2, n_trials = c(10, 0)),
               "at least 1")
  expect_error(sim_blocks("within", n_participants = 2, n_trials = c(10, 7)),
               "must be even")
})

test_that("condition effects can start later in the trial", {
  d <- sim_blocks("within", n_participants = 200, n_trials = 10, n_time = 50,
                  sampling_rate = 25, effect_condition = 1, effect_onset = 1,
                  sd_participant = 0, sd_trial = 0, sd_residual = 0.2,
                  ar1 = 0, seed = 1)
  diff_at <- function(t) {
    x <- d[abs(d$time - t) < 1e-9, ]
    mean(x$y[x$condition == "B"]) - mean(x$y[x$condition == "A"])
  }
  expect_lt(abs(diff_at(0)), 0.05)
  expect_equal(diff_at(1), 0.5, tolerance = 0.1)
  expect_equal(diff_at(1.96), 1, tolerance = 0.05)
  expect_error(sim_blocks(effect_onset = -1), "at least 0")
})

test_that("trial intercepts can be correlated across consecutive trials", {
  d <- sim_blocks("within", "trial", n_participants = 20, n_trials = 400,
                  n_time = 5, sd_participant = 0, sd_slope = 0, sd_trial = 1,
                  sd_residual = 0.01, ar1 = 0, ar1_trials = 0.6, seed = 2)
  d <- d[order(d$participant, d$trial_index), ]
  lag1 <- vapply(split(d$y, d$participant), function(x) {
    stats::cor(x[-1], x[-length(x)])
  }, numeric(1))
  expect_equal(mean(lag1), 0.6, tolerance = 0.03)
  expect_equal(stats::sd(d$y), 1, tolerance = 0.1)
  expect_error(sim_blocks(ar1_trials = 1), "smaller than 1")
})
