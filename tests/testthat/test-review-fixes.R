test_that("data in a different order than the model is refused", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  ds <- d[rev(seq_len(nrow(d))), ]
  expect_error(boot_ci(m, "participant", resample = "wild", B = 9,
                       data = ds), "does not match")
  expect_error(perm_test(m, "condition", "trial", exch_within("participant"),
                         B = 9, data = ds), "does not match")
})

test_that("weights and offsets are refused", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 1)
  d$w <- 1 + seq_len(nrow(d)) %% 3
  m <- lm(y ~ condition, data = d, weights = w)
  expect_error(perm_test(m, "condition", "trial", exch_within("participant"),
                         B = 9), "weights or offsets")
  m2 <- lm(y ~ condition + offset(trial_index), data = d)
  expect_error(boot_ci(m2, "participant", B = 9), "weights or offsets")
})

test_that("series identifiers that repeat across participants are refused", {
  d <- sim_blocks("within", n_participants = 4, n_trials = 6, n_time = 10,
                  seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  expect_error(diag_acf(m, "trial_index", "time"), "more than once")
  expect_error(
    perm_calibrate(m, "condition", "trial", exch_signflip("participant"),
                   ar1 = 0.5, series = "trial_index", time = "time",
                   n_sim = 2, B = 9),
    "more than once"
  )
})

test_that("default bootstrap and LOSO terms exclude the cluster factor", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  bt <- boot_ci(m, "participant", B = 19, seed = 1)
  expect_equal(colnames(bt$boot), "conditionB")
  expect_equal(bt$B_used, 19)
  loso <- diag_loso(m, "participant")
  expect_equal(colnames(loso$estimates), "conditionB")
})

test_that("agreement handles sum-to-zero contrasts", {
  d <- sim_blocks("within", "trial", n_participants = 12, n_trials = 10,
                  effect_condition = 0.6, seed = 2)
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d,
                  contrasts = list(condition = "contr.sum"))
  ag <- diag_agreement(m, "condition", "participant")
  expect_false(ag$flagged)
  expect_equal(ag$estimate, ag$mean_diff, tolerance = 0.05)
})

test_that("groupings that include the cluster are not reported as crossed", {
  d <- sim_blocks("within", n_participants = 6, n_trials = 4, n_time = 5,
                  seed = 3)
  m <- suppressMessages(lme4::lmer(
    y ~ condition + (1 | participant) + (1 | participant:trial_index),
    data = d
  ))
  expect_silent(nested <- nested_grouping(m, d, "participant"))
  expect_length(nested, 0)
})

test_that("the sharp null warns when there is no random slope to remove", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 4)
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d)
  expect_warning(
    perm_calibrate(m, "condition", "trial", exch_signflip("participant"),
                   null = "sharp", n_sim = 2, B = 9, seed = 1),
    "identical"
  )
})

test_that("calibration is reproducible with a non-default RNG kind", {
  skip_on_cran()
  skip_if_not_installed("future")
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 5)
  m <- lm(y ~ condition + participant, data = d)
  c1 <- perm_calibrate(m, "condition", "trial", exch_signflip("participant"),
                       n_sim = 4, B = 9, seed = 2)
  oplan <- future::plan("multisession", workers = 2)
  on.exit(future::plan(oplan), add = TRUE)
  c2 <- perm_calibrate(m, "condition", "trial", exch_signflip("participant"),
                       n_sim = 4, B = 9, seed = 2)
  expect_equal(c1$p_perm, c2$p_perm)
})
