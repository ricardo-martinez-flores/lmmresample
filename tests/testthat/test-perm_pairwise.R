pairwise_data <- function(seed = 1, effect_c = 0.8) {
  d <- sim_blocks("within", "trial", n_participants = 8, n_trials = 12,
                  seed = seed)
  set.seed(seed)
  d$stimulus <- factor(sample(c("A", "B", "C"), nrow(d), replace = TRUE))
  d$y <- d$y + effect_c * (d$stimulus == "C")
  d
}

test_that("perm_pairwise returns one row per pair with the right estimates", {
  d <- pairwise_data()
  m <- lm(y ~ stimulus + participant, data = d)
  res <- perm_pairwise(m, "stimulus", unit = "trial",
                       exchange = exch_signflip("participant"), B = 49,
                       seed = 1)
  expect_s3_class(res, "lmmr_pairwise")
  expect_equal(res$results$comparison, c("B - A", "C - A", "C - B"))
  cf <- stats::coef(m)
  expect_equal(res$results$estimate,
               unname(c(cf["stimulusB"], cf["stimulusC"],
                        cf["stimulusC"] - cf["stimulusB"])),
               tolerance = 1e-8)
  expect_true(all(res$results$p.adjusted >= res$results$p.value - 1e-12))
  expect_equal(res$combine, "max-t")
  expect_s3_class(tidy(res), "data.frame")
  expect_output(print(res), "Pairwise permutation tests")
})

test_that("the reparametrised models reproduce the original fit", {
  d <- pairwise_data()
  m <- lme4::lmer(y ~ stimulus + (1 + stimulus | participant), data = d,
                  control = lme4::lmerControl(calc.derivs = FALSE))
  res <- suppressWarnings(perm_pairwise(
    m, "stimulus", unit = "trial", exchange = exch_signflip("participant"),
    B = 19, seed = 1))
  for (fit in res$models) {
    expect_equal(as.numeric(stats::logLik(fit)),
                 as.numeric(stats::logLik(m)), tolerance = 1e-3)
  }
  cf <- lme4::fixef(m)
  expect_equal(res$results$estimate[1], unname(cf["stimulusB"]),
               tolerance = 1e-4)
  expect_equal(res$results$estimate[3],
               unname(cf["stimulusC"] - cf["stimulusB"]), tolerance = 1e-4)
})

test_that("holm and none adjustments follow their definitions", {
  d <- pairwise_data()
  m <- lm(y ~ stimulus + participant, data = d)
  ex <- exch_signflip("participant")
  holm <- perm_pairwise(m, "stimulus", unit = "trial", exchange = ex,
                        adjust = "holm", B = 49, seed = 3)
  none <- perm_pairwise(m, "stimulus", unit = "trial", exchange = ex,
                        adjust = "none", B = 49, seed = 3)
  expect_equal(holm$results$p.adjusted,
               stats::p.adjust(holm$results$p.value, "holm"))
  expect_equal(none$results$p.adjusted, none$results$p.value)
  expect_equal(holm$results$p.value, none$results$p.value)
})

test_that("a selected subset of pairs is tested", {
  d <- pairwise_data()
  m <- lm(y ~ stimulus + participant, data = d)
  res <- perm_pairwise(m, "stimulus", unit = "trial",
                       exchange = exch_signflip("participant"),
                       pairs = list(c("A", "C"), c("B", "C")), B = 19,
                       seed = 1)
  expect_equal(res$results$comparison, c("C - A", "C - B"))
})

test_that("an effect of another level does not leak into a null pair", {
  d <- pairwise_data(effect_c = 3)
  m <- lm(y ~ stimulus + participant, data = d)
  res <- perm_pairwise(m, "stimulus", unit = "trial",
                       exchange = exch_signflip("participant"), B = 99,
                       seed = 2)
  expect_gt(res$results$p.value[res$results$comparison == "B - A"], 0.05)
  expect_lt(res$results$p.adjusted[res$results$comparison == "C - A"], 0.05)
})

test_that("interactions with time are compared with min-p", {
  d <- sim_blocks("within", n_participants = 8, n_trials = 10, n_time = 8,
                  sampling_rate = 10, prop_condition = 0.5, seed = 4)
  tr <- unique(d["trial"])
  set.seed(4)
  tr$stim <- factor(sample(c("A", "B", "C"), nrow(tr), replace = TRUE))
  d$stim <- tr$stim[match(d$trial, tr$trial)]
  m <- lme4::lmer(y ~ stim * splines::ns(time, 2) + (1 | participant),
                  data = d)
  res <- perm_pairwise(m, "stim", term = "stim:splines::ns(time, 2)",
                       unit = "trial", exchange = exch_signflip("participant"),
                       B = 19, seed = 1)
  expect_equal(res$combine, "min-p")
  expect_true(all(is.na(res$results$estimate)))
  expect_true(all(grepl("chi2", res$results$stat)))
})

test_that("bootstrap intervals are added on request", {
  d <- pairwise_data()
  m <- lm(y ~ stimulus + participant, data = d)
  res <- perm_pairwise(m, "stimulus", unit = "trial",
                       exchange = exch_signflip("participant"), B = 49,
                       seed = 1, cluster = "participant")
  expect_true(all(res$results$conf.low < res$results$estimate))
  expect_true(all(res$results$conf.high > res$results$estimate))
})

test_that("invalid input is refused with clear errors", {
  d <- pairwise_data()
  d$x <- stats::rnorm(nrow(d))
  m <- lm(y ~ stimulus + x + participant, data = d)
  ex <- exch_signflip("participant")
  expect_error(perm_pairwise(m, "stimulus", term = "x", unit = "trial",
                             exchange = ex, B = 9),
               "does not contain")
  expect_error(perm_pairwise(m, "x", unit = "trial", exchange = ex, B = 9),
               "factor or a character")
  expect_error(perm_pairwise(m, "stimulus", unit = "trial", exchange = ex,
                             pairs = list(c("A", "Z")), B = 9),
               "Unknown level")
  expect_error(perm_pairwise(m, "stimulus", unit = "trial", exchange = ex,
                             pairs = list(c("A", "A")), B = 9),
               "two different")
})
