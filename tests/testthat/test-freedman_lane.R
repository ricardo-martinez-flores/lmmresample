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
  eng <- build_engine(m, d, "time", test, "trial", exch_within("participant"),
                      time = "time")
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
              method = "freedman-lane", time = "time", B = 9),
    "equal size"
  )
  expect_error(
    perm_test(m, "time", "trial", exch_within("participant"),
              method = "freedman-lane", B = 9),
    "requires `time`"
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

test_that("random effects at or below the units do not absorb the effect", {
  # Between-group effect with participants as the permuted units: the
  # participant intercepts must not be kept in the reduced fitted values.
  set.seed(1)
  n <- 24
  part <- data.frame(participant = factor(seq_len(n)),
                     group = factor(rep(c("a", "b", "c"), length.out = n)))
  part$u <- stats::rnorm(n, 0, 0.5)
  d <- part[rep(seq_len(n), each = 10), ]
  d$x <- stats::rnorm(nrow(d))
  d$idx <- rep(1:10, n)
  d$y <- c(a = 0, b = 1, c = 2)[as.character(d$group)] + d$u + 0.3 * d$x +
    stats::rnorm(nrow(d))
  m <- lme4::lmer(y ~ group + x + (1 | participant), data = d)
  res <- perm_test(m, "group", "participant", exch_free(), time = "idx",
                   B = 99, seed = 1)
  expect_equal(res$method, "freedman-lane")
  expect_lt(res$p.value, 0.05)
  expect_lt(stats::median(res$null, na.rm = TRUE), 10)

  # Trial intercepts below the participant blocks
  d2 <- sim_blocks("within", n_participants = 10, n_trials = 12, n_time = 5,
                   seed = 2)
  tr <- levels(d2$trial)
  c3 <- stats::setNames(rep(c("a", "b", "c"), length.out = length(tr)), tr)
  d2$cond3 <- factor(c3[as.character(d2$trial)])
  d2$y <- d2$y + c(a = 0, b = 0.8, c = 1.6)[as.character(d2$cond3)]
  m2 <- suppressMessages(lme4::lmer(y ~ cond3 + (1 | participant) +
                                      (1 | trial), data = d2))
  res2 <- suppressWarnings(perm_test(m2, "cond3", "trial",
                                     exch_within("participant"),
                                     time = "time", B = 49, seed = 1))
  expect_lt(res2$p.value, 0.05)
})

test_that("residuals are matched by time, not by row order", {
  d <- fl_data(ar1 = 0.9)
  set.seed(4)
  ds <- d[sample(nrow(d)), ]
  rownames(ds) <- NULL
  m <- lm(y ~ condition + time + participant, data = ds)
  test <- resolve_test(m, ds, "time", NULL, "freedman-lane")
  eng <- build_engine(m, ds, "time", test, "trial",
                      exch_within("participant"), time = "time")
  set.seed(1)
  draws <- eng$draw(3)
  fit_r <- stats::lm(y ~ condition + participant, data = ds)
  e <- unname(ds$y - stats::fitted(fit_r))
  e_new <- unname(eng$make_data(draws, 1)$y - stats::fitted(fit_r))
  series <- function(v) {
    lapply(split(seq_along(v), ds$trial), function(i) v[i][order(ds$time[i])])
  }
  orig <- series(e)
  new <- series(e_new)
  matched <- vapply(new, function(x) {
    any(vapply(orig, function(y) isTRUE(all.equal(x, y)), logical(1)))
  }, logical(1))
  expect_true(all(matched))
})

test_that("random-effects terms left empty are dropped", {
  d <- fl_data()
  m <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ condition + time + (1 | participant) +
                 (0 + condition | participant), data = d)
  ))
  f <- paste(deparse(reduced_formula(m, "condition")), collapse = "")
  expect_false(grepl("(0 |", f, fixed = TRUE))
  expect_match(f, "(1 | participant)", fixed = TRUE)
})
