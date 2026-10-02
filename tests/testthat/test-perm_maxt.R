maxt_data <- function(seed = 1, effect = 0, ...) {
  d <- sim_blocks("within", "trial", n_participants = 8, n_trials = 8,
                  n_time = 10, effect_condition = effect, seed = seed, ...)
  set.seed(seed)
  d$noise <- stats::rnorm(nrow(d))
  d
}

test_that("perm_maxt returns a complete lmmr_maxt object", {
  d <- maxt_data(effect = 1)
  res <- perm_maxt(
    mean = perm_spec(lm(y ~ condition + participant, data = d), "condition"),
    peak = perm_spec(lm(peak ~ condition + participant, data = d), "condition"),
    noise = perm_spec(lm(noise ~ condition + participant, data = d),
                      "condition"),
    unit = "trial", exchange = exch_within("participant"), B = 99, seed = 1
  )
  expect_s3_class(res, "lmmr_maxt")
  expect_equal(res$results$test, c("mean", "peak", "noise"))
  expect_equal(dim(res$null), c(99, 3))
  expect_true(all(res$results$p.adjusted >= res$results$p.value))
  expect_lt(res$results$p.adjusted[1], 0.05)
  expect_gt(res$results$p.adjusted[3], 0.05)
})

test_that("each column of the joint null matches a separate perm_test", {
  d <- maxt_data()
  m1 <- lm(y ~ condition + participant, data = d)
  m2 <- lm(noise ~ condition + participant, data = d)
  ex <- exch_within("participant")
  joint <- perm_maxt(a = perm_spec(m1, "condition"),
                     b = perm_spec(m2, "condition"),
                     unit = "trial", exchange = ex, B = 30, seed = 5)
  single <- perm_test(m1, "condition", "trial", ex, B = 30, seed = 5)
  expect_equal(joint$null[, "a"], abs(single$null))
  expect_equal(joint$results$p.value[1], single$p.value)
})

test_that("max-t adjustment follows the definition", {
  null <- cbind(c(1, 2, 3, 4), c(4, 1, 1, 1))
  observed <- c(3.5, 1.5)
  ss <- maxt_adjust(observed, null, "single-step")
  expect_equal(ss, c((2 + 1) / 5, (4 + 1) / 5))
  sd <- maxt_adjust(observed, null, "step-down")
  # first (largest) uses the maximum over both tests; second only its own
  expect_equal(sd[1], (2 + 1) / 5)
  expect_equal(sd[2], max(sd[1], (1 + 1) / 5))
  expect_true(all(sd <= ss))
})

test_that("step-down p-values are monotone in the observed statistics", {
  set.seed(3)
  null <- matrix(abs(stats::rnorm(400)), ncol = 4)
  observed <- c(2.5, 0.5, 1.8, 3)
  p <- maxt_adjust(observed, null, "step-down")
  ord <- order(observed, decreasing = TRUE)
  expect_true(all(diff(p[ord]) >= 0))
})

test_that("the same relabelling is shared across models with different data", {
  d <- maxt_data()
  d2 <- d[sample(nrow(d)), ]
  m1 <- lm(y ~ condition + participant, data = d)
  m2 <- lm(noise ~ condition + participant, data = d2)
  res <- perm_maxt(a = perm_spec(m1, "condition"),
                   b = perm_spec(m2, "condition"),
                   unit = "trial", exchange = exch_signflip("participant"),
                   B = 20, seed = 2)
  expect_equal(res$B_used, 20)
})

test_that("invalid families are refused", {
  d <- maxt_data()
  ex <- exch_within("participant")
  m <- lm(y ~ condition + participant, data = d)
  s <- perm_spec(m, "condition")
  expect_error(perm_maxt(a = s, unit = "trial", exchange = ex), "at least two")
  expect_error(perm_maxt(s, s, unit = "trial", exchange = ex), "unique names")
  expect_error(perm_maxt(a = s, b = m, unit = "trial", exchange = ex),
               "perm_spec")
  m_p <- lm(y ~ condition + participant, data = d)
  expect_error(
    perm_maxt(a = s, b = perm_spec(m_p, "participant", coef = "participantp2"),
              unit = "trial", exchange = ex),
    "same variable"
  )
  d_sub <- d[d$participant != levels(d$participant)[1], ]
  m_sub <- lm(y ~ condition, data = d_sub)
  expect_error(
    perm_maxt(a = s, b = perm_spec(m_sub, "condition"), unit = "trial",
              exchange = ex),
    "same units"
  )
})

test_that("print, tidy and plot methods work", {
  d <- maxt_data()
  res <- perm_maxt(
    a = perm_spec(lm(y ~ condition + participant, data = d), "condition"),
    b = perm_spec(lm(noise ~ condition + participant, data = d), "condition"),
    unit = "trial", exchange = exch_within("participant"), B = 19, seed = 1
  )
  expect_output(print(res), "max-t")
  expect_output(print(res), "family-wise")
  expect_equal(nrow(tidy(res)), 2)
  expect_equal(nrow(tidy(res, type = "null")), 2 * 19)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(res), "ggplot")
})

test_that("max-t controls the family-wise error rate", {
  skip_on_cran()
  n_sim <- 60
  any_unadj <- any_adj <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    d <- sim_blocks("within", "trial", n_participants = 10, n_trials = 8,
                    seed = i)
    set.seed(i)
    for (j in 1:4) d[[paste0("f", j)]] <- d$y + stats::rnorm(nrow(d))
    specs <- lapply(1:4, function(j) {
      perm_spec(lm(stats::as.formula(paste0("f", j,
                                            " ~ condition + participant")),
                   data = d), "condition")
    })
    names(specs) <- paste0("f", 1:4)
    res <- do.call(perm_maxt, c(specs, list(unit = "trial",
                                            exchange = exch_signflip("participant"),
                                            B = 99, seed = i)))
    any_unadj[i] <- any(res$results$p.value <= 0.05)
    any_adj[i] <- any(res$results$p.adjusted <= 0.05)
  }
  expect_lt(mean(any_adj), 0.15)
  expect_gte(mean(any_unadj), mean(any_adj))
})

test_that("perm_maxt supports Freedman-Lane with interactions", {
  d <- sim_blocks("within", n_participants = 8, n_trials = 6, n_time = 10,
                  effect_condition = 1, effect_onset = 0.1, sd_residual = 0.3,
                  seed = 1)
  set.seed(1)
  d$noise <- stats::rnorm(nrow(d))
  m1 <- lme4::lmer(y ~ condition * time + (1 | participant), data = d)
  m2 <- lme4::lmer(noise ~ condition * time + (1 | participant), data = d)
  res <- perm_maxt(signal = perm_spec(m1, "condition:time"),
                   noise = perm_spec(m2, "condition:time"),
                   unit = "trial", exchange = exch_signflip("participant"),
                   B = 49, seed = 2)
  expect_equal(res$method, "freedman-lane")
  expect_lt(res$results$p.adjusted[1], 0.05)
  expect_gt(res$results$p.adjusted[2], 0.05)
  expect_output(print(res), "Freedman-Lane")

  single <- perm_test(m1, "condition:time", "trial",
                      exch_signflip("participant"), B = 49, seed = 2)
  expect_equal(res$null[, "signal"], abs(single$null))
})

test_that("Freedman-Lane families can mix terms and data orders", {
  d <- sim_blocks("within", n_participants = 6, n_trials = 4, n_time = 8,
                  seed = 3)
  d2 <- d[sample(nrow(d)), ]
  m1 <- suppressMessages(
    lme4::lmer(y ~ condition + time + (1 | participant), data = d)
  )
  m2 <- suppressMessages(
    lme4::lmer(y ~ condition * time + (1 | participant), data = d2)
  )
  res <- perm_maxt(main = perm_spec(m1, "condition"),
                   inter = perm_spec(m2, "condition:time"),
                   unit = "trial", exchange = exch_signflip("participant"),
                   B = 19, seed = 1)
  expect_equal(res$method, "freedman-lane")
  expect_equal(res$B_used, 19)

  m3 <- suppressMessages(
    lme4::lmer(y ~ condition * splines::ns(time, 2) + (1 | participant),
               data = d)
  )
  expect_error(
    perm_maxt(a = perm_spec(m1, "condition"),
              b = perm_spec(m3, "condition:splines::ns(time, 2)"),
              unit = "trial", exchange = exch_signflip("participant"), B = 9),
    "single coefficient"
  )
  expect_error(
    perm_maxt(a = perm_spec(m1, "condition"), b = perm_spec(m1, "time"),
              unit = "trial", exchange = exch_signflip("participant"),
              method = "relabel", B = 9),
    "same variable"
  )
})
