fl_data <- function(seed = 1, ...) {
  sim_blocks("within", "timeseries", n_participants = 8, n_trials = 8,
             n_time = 10, seed = seed, ...)
}

test_that("the reduced formula drops the term and its random slopes", {
  d <- fl_data()
  m <- lme4::lmer(y ~ condition * time + (1 + condition | participant),
                  data = d)
  f <- reduced_formula(m, "condition:time")
  expect_false("condition:time" %in% attr(stats::terms(lme4::nobars(f)),
                                           "term.labels"))
  expect_match(paste(deparse(f), collapse = ""), "condition | participant",
               fixed = TRUE)

  m2 <- lme4::lmer(y ~ condition + time + (1 + condition | participant),
                   data = d)
  f2 <- paste(deparse(reduced_formula(m2, "condition")), collapse = "")
  expect_match(f2, "(1 | participant)", fixed = TRUE)
  expect_false(grepl("condition", f2))

  m3 <- lm(y ~ condition + time + participant, data = d)
  f3 <- reduced_formula(m3, "time")
  expect_equal(attr(stats::terms(f3), "term.labels"),
               c("condition", "participant"))
})

test_that("method selection follows the term", {
  d <- fl_data()
  m <- lme4::lmer(y ~ condition * time + (1 | participant), data = d)
  expect_equal(resolve_test(m, d, "condition:time", NULL, "auto")$method,
               "freedman-lane")
  m2 <- lme4::lmer(y ~ condition + time + (1 | participant), data = d)
  expect_equal(resolve_test(m2, d, "condition", NULL, "auto")$method,
               "relabel")
  expect_error(resolve_test(m, d, "condition", NULL, "auto"),
               "higher-order term")
  d$high <- as.integer(d$y > 0)
  g <- glm(high ~ condition + time, family = binomial, data = d)
  expect_error(resolve_test(g, d, "time", NULL, "freedman-lane"), "Gaussian")
})

test_that("Freedman-Lane tests an interaction with sign-flipping", {
  d <- fl_data(effect_condition = 1, effect_onset = 0.1, sd_residual = 0.3)
  m <- lme4::lmer(y ~ condition * time + (1 | participant), data = d)
  res <- perm_test(m, "condition:time", "trial",
                   exch_signflip("participant"), B = 49, seed = 1)
  expect_equal(res$method, "freedman-lane")
  expect_equal(res$coef, "conditionB:time")
  expect_lt(res$p.value, 0.05)
  expect_output(print(res), "Freedman-Lane")
})

test_that("terms with several coefficients use a joint chi-square", {
  d <- fl_data()
  m <- lme4::lmer(y ~ condition * splines::ns(time, 3) + (1 | participant),
                  data = d)
  res <- perm_test(m, "condition:splines::ns(time, 3)", "trial",
                   exch_signflip("participant"), B = 19, seed = 1)
  expect_equal(res$stat_type, "chi2")
  expect_equal(res$df, 3)
  expect_true(all(res$null[!is.na(res$null)] >= 0))
  expect_output(print(res), "chi2\\(3\\)")
  expect_error(
    perm_test(m, "condition:splines::ns(time, 3)", "trial",
              exch_signflip("participant"), alternative = "greater", B = 9),
    "One-sided"
  )
})

test_that("residual permutation between units keeps whole trajectories", {
  d <- fl_data()
  m <- lm(y ~ condition + time + participant, data = d)
  test <- resolve_test(m, d, "time", NULL, "freedman-lane")
  eng <- build_engine(m, d, "time", test, "trial", exch_within("participant"))
  set.seed(1)
  draws <- eng$draw(5)
  fit_r <- stats::lm(y ~ condition + participant, data = d)
  e <- unname(d$y - stats::fitted(fit_r))
  new <- eng$make_data(draws, 1)
  e_new <- unname(new$y - stats::fitted(fit_r))
  by_trial <- split(e, d$trial)
  by_trial_new <- split(e_new, d$trial)
  matched <- vapply(by_trial_new, function(x) {
    any(vapply(by_trial, function(y) isTRUE(all.equal(x, y)), logical(1)))
  }, logical(1))
  expect_true(all(matched))
  expect_equal(sort(e_new), sort(e))
})

test_that("residual permutation requires units of equal size", {
  d <- fl_data()
  d <- d[!(d$trial == levels(d$trial)[1] & d$time > 0.1), ]
  m <- lm(y ~ condition + time + participant, data = d)
  expect_error(
    perm_test(m, "time", "trial", exch_within("participant"),
              method = "freedman-lane", B = 9),
    "equal size"
  )
  expect_s3_class(
    perm_test(m, "time", "trial", exch_signflip("participant"),
              method = "freedman-lane", B = 9, seed = 1),
    "lmmr_perm"
  )
})

test_that("Freedman-Lane controls type I error for an interaction", {
  skip_on_cran()
  n_sim <- 40
  rej <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    d <- sim_blocks("within", n_participants = 10, n_trials = 8, n_time = 15,
                    effect_condition = 0.5, sd_slope = 0.4, ar1 = 0.9,
                    seed = i)
    m <- lme4::lmer(y ~ condition * time + (1 | participant), data = d)
    rej[i] <- perm_test(m, "condition:time", "trial",
                        exch_signflip("participant"), B = 49,
                        seed = i)$p.value <= 0.05
  }
  expect_lt(mean(rej), 0.15)
})

test_that("calibration works with Freedman-Lane", {
  d <- fl_data()
  m <- lme4::lmer(y ~ condition * time + (1 | participant), data = d)
  cal <- perm_calibrate(m, "condition:time", "trial",
                        exch_signflip("participant"), n_sim = 4, B = 9,
                        seed = 1)
  expect_s3_class(cal, "lmmr_calib")
  expect_equal(cal$method, "freedman-lane")
  expect_length(cal$p_perm, 4)
})
