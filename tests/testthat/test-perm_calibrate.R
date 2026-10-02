calib_data <- function(seed = 1, ...) {
  sim_blocks("within", "trial", n_participants = 8, n_trials = 8,
             n_time = 10, seed = seed, ...)
}

test_that("perm_calibrate returns a complete lmmr_calib object", {
  d <- calib_data()
  m <- lm(y ~ condition + participant, data = d)
  cal <- perm_calibrate(m, "condition", "trial", exch_within("participant"),
                        n_sim = 10, B = 19, seed = 1)

  expect_s3_class(cal, "lmmr_calib")
  expect_length(cal$p_perm, 10)
  expect_length(cal$p_wald, 10)
  expect_equal(nrow(cal$rejection), 2)
  expect_true(all(cal$p_perm > 0 & cal$p_perm <= 1))
  expect_true(all(cal$B_used == 19))
})

test_that("calibration is reproducible with a seed", {
  d <- calib_data()
  m <- lm(y ~ condition + participant, data = d)
  c1 <- perm_calibrate(m, "condition", "trial", exch_within("participant"),
                       n_sim = 5, B = 9, seed = 2)
  c2 <- perm_calibrate(m, "condition", "trial", exch_within("participant"),
                       n_sim = 5, B = 9, seed = 2)
  expect_identical(c1$p_perm, c2$p_perm)
})

test_that("the null simulator removes the tested effect", {
  d <- calib_data(effect_condition = 2)
  m <- lm(y ~ condition + participant, data = d)
  sim <- make_null_simulator(m, d, "conditionB", "mean")
  set.seed(1)
  diffs <- replicate(100, {
    y <- sim()
    mean(y[d$condition == "B"]) - mean(y[d$condition == "A"])
  })
  expect_lt(abs(mean(diffs)), 0.1)
})

test_that("the sharp null removes the random slope of the tested term", {
  d <- sim_blocks("within", "trial", n_participants = 20, n_trials = 10,
                  sd_slope = 1, seed = 2)
  m <- suppressMessages(
    lme4::lmer(y ~ condition + (1 + condition | participant), data = d)
  )
  theta <- lme4::getME(m, "theta")
  zeroed <- zero_theta_rows(m, theta, "conditionB")
  expect_equal(zeroed[[1]], theta[[1]])
  expect_true(all(zeroed[-1] == 0))

  spread <- function(null) {
    sim <- make_null_simulator(m, d, "conditionB", null)
    set.seed(3)
    stats::sd(replicate(100, {
      y <- sim()
      mean(y[d$condition == "B"]) - mean(y[d$condition == "A"])
    }))
  }
  expect_gt(spread("mean"), 1.5 * spread("sharp"))
})

test_that("AR(1) residuals have the requested autocorrelation", {
  d <- sim_blocks("within", "timeseries", n_participants = 2, n_trials = 20,
                  n_time = 500, ar1 = 0, seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  sim <- make_null_simulator(m, d, "conditionB", "mean", ar1 = 0.8,
                             series = "trial", time = "time")
  set.seed(4)
  y <- sim()
  r <- y - stats::fitted(m) + stats::coef(m)["conditionB"] *
    (d$condition == "B")
  lag1 <- vapply(split(r, d$trial), function(x) {
    stats::cor(x[-1], x[-length(x)])
  }, numeric(1))
  expect_equal(mean(lag1), 0.8, tolerance = 0.03)
})

test_that("invalid calibrations are refused", {
  d <- calib_data()
  ex <- exch_within("participant")
  m <- lm(y ~ condition + participant, data = d)

  expect_error(perm_calibrate(m, "condition", "trial", ex, ar1 = 0.5),
               "series")
  expect_error(perm_calibrate(m, "condition", "trial", ex, ar1 = 1),
               "smaller than 1")
  expect_error(perm_calibrate(m, "condition", "trial", ex, alpha = 0),
               "greater than 0")

  m_log <- lm(log(y + 10) ~ condition, data = d)
  expect_error(perm_calibrate(m_log, "condition", "trial", ex),
               "column of the data")

})

test_that("print, tidy and plot methods work", {
  d <- calib_data()
  m <- lm(y ~ condition + participant, data = d)
  cal <- perm_calibrate(m, "condition", "trial", exch_within("participant"),
                        n_sim = 10, B = 19, seed = 1)
  expect_output(print(cal), "Rejection rate at alpha")
  expect_equal(nrow(tidy(cal)), 2)
  expect_equal(nrow(tidy(cal, type = "pvalues")), 10)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(cal), "ggplot")
  expect_s3_class(plot(cal, type = "rejection"), "ggplot")
})

test_that("calibration detects the liberal Wald test with autocorrelation", {
  skip_on_cran()
  d <- sim_blocks("within", "timeseries", n_participants = 6, n_trials = 8,
                  n_time = 30, ar1 = 0.95, seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  cal <- perm_calibrate(m, "condition", "trial", exch_within("participant"),
                        ar1 = 0.95, series = "trial", time = "time",
                        n_sim = 60, B = 49, seed = 1)
  rej <- tidy(cal)
  expect_gt(rej$rejection[rej$method == "Wald"], 0.2)
  expect_lt(rej$rejection[rej$method == "Permutation"], 0.15)
})
