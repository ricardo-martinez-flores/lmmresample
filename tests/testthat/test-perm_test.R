within_trials <- function(seed = 1, effect = 0, ...) {
  sim_blocks("within", "trial", n_participants = 8, n_trials = 8,
             n_time = 10, effect_condition = effect, seed = seed, ...)
}

test_that("perm_test returns a complete lmmr_perm object", {
  d <- within_trials(effect = 1)
  m <- lm(y ~ condition + participant, data = d)
  res <- perm_test(m, "condition", "trial", exch_within("participant"),
                   B = 99, seed = 1)

  expect_s3_class(res, "lmmr_perm")
  expect_equal(res$coef, "conditionB")
  expect_equal(res$statistic,
               unname(summary(m)$coefficients["conditionB", "t value"]))
  expect_length(res$null, 99)
  expect_equal(res$B_used, 99)
  expect_true(res$p.value > 0 && res$p.value <= 1)
  expect_equal(res$p.value * 100, round(res$p.value * 100))
})

test_that("a strong effect gives the smallest attainable p-value", {
  d <- within_trials(effect = 3, sd_trial = 0.1)
  m <- lm(y ~ condition + participant, data = d)
  res <- perm_test(m, "condition", "trial", exch_within("participant"),
                   B = 99, seed = 1)
  expect_equal(res$p.value, 1 / 100)
})

test_that("results are reproducible with a seed", {
  d <- within_trials()
  m <- lm(y ~ condition + participant, data = d)
  r1 <- perm_test(m, "condition", "trial", exch_within("participant"),
                  B = 50, seed = 3)
  r2 <- perm_test(m, "condition", "trial", exch_within("participant"),
                  B = 50, seed = 3)
  r3 <- perm_test(m, "condition", "trial", exch_within("participant"),
                  B = 50, seed = 4)
  expect_identical(r1$null, r2$null)
  expect_false(identical(r1$null, r3$null))
})

test_that("results do not depend on the parallel backend", {
  skip_on_cran()
  skip_if_not_installed("future")
  d <- within_trials()
  m <- lm(y ~ condition + participant, data = d)
  seq_res <- perm_test(m, "condition", "trial", exch_within("participant"),
                       B = 40, seed = 7)
  old <- future::plan("multisession", workers = 2)
  on.exit(future::plan(old), add = TRUE)
  par_res <- perm_test(m, "condition", "trial", exch_within("participant"),
                       B = 40, seed = 7)
  expect_equal(par_res$null, seq_res$null)
})

test_that("one-sided alternatives use the correct tail", {
  null <- c(-2, -1, 0, 1, 2)
  expect_equal(perm_pvalue(1.5, null, "greater"), 2 / 6)
  expect_equal(perm_pvalue(1.5, null, "less"), 5 / 6)
  expect_equal(perm_pvalue(1.5, null, "two.sided"), 3 / 6)
  expect_equal(perm_pvalue(-1.5, null, "two.sided"), 3 / 6)
  expect_true(is.na(perm_pvalue(1, numeric(0), "two.sided")))
})

test_that("lm and lmer models are supported; generalized models are not", {
  d <- within_trials(effect = 0.5)
  d$high <- as.integer(d$y > stats::median(d$y))
  ex <- exch_within("participant")

  r_lm <- perm_test(lm(y ~ condition + participant, data = d),
                    "condition", "trial", ex, B = 19, seed = 1)
  r_lmer <- perm_test(lme4::lmer(y ~ condition + (1 | participant), data = d),
                      "condition", "trial", ex, B = 19, seed = 1)
  expect_equal(r_lm$stat_label, "t")
  expect_equal(r_lmer$stat_label, "t")
  expect_true(r_lm$B_used > 0 && r_lmer$B_used > 0)

  expect_error(perm_test(glm(high ~ condition, family = binomial, data = d),
                         "condition", "trial", ex), "linear model")
  g <- suppressMessages(lme4::glmer(high ~ condition + (1 | participant),
                                    family = binomial, data = d))
  expect_error(perm_test(g, "condition", "trial", ex), "linear model")
})

test_that("between-group designs permute participants", {
  d <- sim_blocks("between", "trial", n_participants = 10, n_trials = 6,
                  effect_group = 0.5, seed = 2)
  m <- lme4::lmer(y ~ group + (1 | participant), data = d)
  res <- perm_test(m, "group", "participant", exch_free(), B = 49, seed = 1)
  expect_equal(res$n_units, 10)
  expect_equal(res$coef, "grouptreatment")
})

test_that("invalid tests are refused with informative errors", {
  d <- within_trials()
  d$three <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  ex <- exch_within("participant")

  m <- lm(y ~ condition + participant, data = d)
  expect_error(perm_test(m, "group", "trial", ex), "not a fixed-effect term")
  expect_error(perm_test(m, "condition", "participant", exch_free()),
               "constant within each")
  expect_error(perm_test(m, "condition", "trial", ex, B = 0), "at least 1")
  expect_error(perm_test(m, "condition", "trial", ex, coef = "x"),
               "not a coefficient")

  m_int <- lm(y ~ condition * trial_index, data = d)
  expect_error(perm_test(m_int, "condition", "trial", ex), "interaction")


  expect_error(perm_test(list(), "condition", "trial", ex), "linear model")

  m_sub <- lm(y ~ condition, data = d, subset = trial_index > 2)
  expect_error(perm_test(m_sub, "condition", "trial", ex), "rows")
})

test_that("a factor with several levels can be tested through one coefficient", {
  d <- within_trials()
  d$three <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  m <- lm(y ~ three, data = d)
  res <- perm_test(m, "three", "trial", exch_free(), coef = "threec",
                   B = 19, seed = 1)
  expect_equal(res$coef, "threec")
})

test_that("data that cannot be recovered can be supplied", {
  d <- within_trials()
  m <- lm(y ~ condition + participant, data = d)
  m$call$data <- quote(object_that_does_not_exist)
  expect_error(
    perm_test(m, "condition", "trial", exch_within("participant"), B = 19),
    "Cannot recover the data"
  )
  res <- perm_test(m, "condition", "trial", exch_within("participant"),
                   B = 19, seed = 1, data = d)
  expect_s3_class(res, "lmmr_perm")
})

test_that("small permutation spaces trigger a warning", {
  d <- sim_blocks("between", "participant", n_participants = 6, seed = 1)
  m <- lm(y ~ group, data = d)
  expect_warning(
    perm_test(m, "group", "participant", exch_free(), B = 99, seed = 1),
    "distinct relabellings"
  )
})

test_that("blocks without variation trigger a warning", {
  d <- within_trials()
  d$condition[d$participant == levels(d$participant)[1]] <- "A"
  m <- lm(y ~ condition + participant, data = d)
  expect_warning(
    perm_test(m, "condition", "trial", exch_within("participant"), B = 19,
              seed = 1),
    "do not contribute"
  )
})

test_that("print, summary and tidy methods work", {
  d <- within_trials()
  m <- lm(y ~ condition + participant, data = d)
  res <- perm_test(m, "condition", "trial", exch_within("participant"),
                   B = 19, seed = 1)

  expect_output(print(res), "Permutation test for `condition`")
  expect_output(print(res), "exploratory")
  expect_output(print(summary(res)), "Refit status")

  td <- tidy(res)
  expect_equal(nrow(td), 1)
  expect_named(td, c("term", "coef", "statistic", "df", "p.value",
                     "alternative", "B", "B_used", "method"))
  expect_equal(nrow(tidy(res, type = "null")), 19)
})

test_that("plots are ggplot objects", {
  skip_if_not_installed("ggplot2")
  d <- within_trials()
  m <- lm(y ~ condition + participant, data = d)
  res <- perm_test(m, "condition", "trial", exch_within("participant"),
                   B = 19, seed = 1)
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, type = "trace"), "ggplot")
})

test_that("type I error is controlled with autocorrelated time series", {
  skip_on_cran()
  # Under the null, a naive Wald test on autocorrelated samples is liberal,
  # whereas permuting whole trials keeps the nominal rate.
  n_sim <- 40
  wald <- perm <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    d <- sim_blocks("within", "timeseries", n_participants = 6, n_trials = 8,
                    n_time = 30, ar1 = 0.95, sd_trial = 0.3, seed = i)
    m <- lm(y ~ condition + participant, data = d)
    wald[i] <- summary(m)$coefficients["conditionB", 4] < 0.05
    perm[i] <- perm_test(m, "condition", "trial",
                         exch_within("participant"), B = 99,
                         seed = i)$p.value < 0.05
  }
  expect_gt(mean(wald), 0.3)
  expect_lt(mean(perm), 0.15)
})

test_that("perm_test supports sign-flipping within participants", {
  d <- within_trials(effect = 0.5)
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d)
  res <- perm_test(m, "condition", "trial", exch_signflip("participant"),
                   B = 49, seed = 1)
  expect_s3_class(res, "lmmr_perm")
  expect_output(print(res), "sign-flipped")

  d$condition[d$participant == levels(d$participant)[1]] <- "A"
  m2 <- lm(y ~ condition + participant, data = d)
  expect_error(
    perm_test(m2, "condition", "trial", exch_signflip("participant"), B = 9),
    "both values in every block"
  )
})

test_that("fast refits give the same statistics as default refits", {
  d <- within_trials(effect = 0.5)
  m <- lme4::lmer(y ~ condition + (1 + condition | participant), data = d)
  fast <- suppressWarnings(perm_test(m, "condition", "trial",
                                     exch_signflip("participant"), B = 9,
                                     seed = 1))
  old <- options(lmmresample.fast = FALSE)
  on.exit(options(old), add = TRUE)
  slow <- suppressWarnings(perm_test(m, "condition", "trial",
                                     exch_signflip("participant"), B = 9,
                                     seed = 1))
  ok <- !is.na(fast$null) & !is.na(slow$null)
  expect_equal(fast$null[ok], slow$null[ok], tolerance = 1e-3)
})
