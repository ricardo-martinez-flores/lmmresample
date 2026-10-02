test_that("diag_acf recovers the residual autocorrelation", {
  d <- sim_blocks("within", n_participants = 6, n_trials = 30, n_time = 200,
                  ar1 = 0.8, sd_trial = 0, seed = 1)
  m <- lm(y ~ condition + participant + factor(time), data = d)
  ac <- diag_acf(m, series = "trial", time = "time")
  expect_s3_class(ac, "lmmr_acf")
  expect_equal(ac$phi, 0.8, tolerance = 0.05)
  expect_equal(nrow(tidy(ac)), 20)
  expect_true(all(diff(ac$acf$acf[1:5]) < 0))
  expect_output(print(ac), "strongly autocorrelated")
})

test_that("diag_acf output can be passed to perm_calibrate", {
  d <- sim_blocks("within", n_participants = 6, n_trials = 6, n_time = 20,
                  seed = 2)
  m <- lm(y ~ condition + participant, data = d)
  ac <- diag_acf(m, "trial", "time")
  cal <- perm_calibrate(m, "condition", "trial",
                        exch_signflip("participant"), ar1 = ac,
                        series = "trial", time = "time", n_sim = 3, B = 9,
                        seed = 1)
  expect_equal(cal$ar1, ac$phi)
})

test_that("diag_acf validates its input", {
  d <- sim_blocks("within", n_participants = 4, n_trials = 4, n_time = 5,
                  seed = 1)
  m <- lm(y ~ condition, data = d)
  expect_error(diag_acf(m, "trial", "clock"), "not found")
  expect_error(diag_acf(m, "trial", "time", lag_max = 10), "at least 12")
})

test_that("diag_timecourse detects a missing condition by time interaction", {
  d <- sim_blocks("within", n_participants = 10, n_trials = 10, n_time = 30,
                  effect_condition = 1, effect_onset = 0.3, sd_residual = 0.3,
                  seed = 3)
  m_bad <- lm(y ~ condition + time + participant, data = d)
  m_good <- lm(y ~ condition * factor(time) + participant, data = d)
  tc_bad <- diag_timecourse(m_bad, "time", by = "condition",
                            cluster = "participant")
  tc_good <- diag_timecourse(m_good, "time", by = "condition",
                             cluster = "participant")
  expect_s3_class(tc_bad, "lmmr_timecourse")
  expect_gt(max(tc_bad$max_abs), 0.3)
  expect_lt(max(tc_good$max_abs), 0.1)
  expect_output(print(tc_bad), "flexible")
  expect_equal(sort(unique(tidy(tc_bad)$group)), c("A", "B"))
})

test_that("diag_timecourse bins continuous time", {
  d <- sim_blocks("within", n_participants = 4, n_trials = 4, n_time = 120,
                  sampling_rate = 100, seed = 4)
  m <- lm(y ~ condition + participant, data = d)
  tc <- diag_timecourse(m, "time", bins = 12)
  expect_equal(length(unique(tidy(tc)$time)), 12)
  tc_raw <- diag_timecourse(m, "time", bins = NULL)
  expect_equal(length(unique(tidy(tc_raw)$time)), 120)
})

test_that("diagnostic plots are ggplot objects", {
  skip_if_not_installed("ggplot2")
  d <- sim_blocks("within", n_participants = 4, n_trials = 6, n_time = 30,
                  seed = 5)
  m <- lm(y ~ condition + participant, data = d)
  expect_s3_class(plot(diag_acf(m, "trial", "time")), "ggplot")
  expect_s3_class(plot(diag_timecourse(m, "time", by = "condition")), "ggplot")
  expect_s3_class(plot(diag_timecourse(m, "time")), "ggplot")
})
