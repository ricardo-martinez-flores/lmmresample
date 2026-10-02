test_that("diag_loso matches manual leave-one-out refits", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 8,
                  effect_condition = 0.5, seed = 1)
  m <- lm(y ~ condition + participant, data = d)
  loso <- diag_loso(m, "participant", terms = "conditionB")
  expect_s3_class(loso, "lmmr_loso")
  expect_equal(dim(loso$estimates), c(6, 1))
  id <- levels(d$participant)[3]
  manual <- stats::coef(lm(y ~ condition + participant,
                           data = droplevels(d[d$participant != id, ])))
  expect_equal(unname(loso$estimates[id, "conditionB"]),
               unname(manual["conditionB"]))
  expect_equal(nrow(tidy(loso)), 6)
  expect_output(print(loso), "Most influential")
})

test_that("diag_loso flags an influential participant", {
  d <- sim_blocks("within", "trial", n_participants = 10, n_trials = 10,
                  effect_condition = 0, sd_trial = 0.1, sd_residual = 0.1,
                  seed = 2)
  out <- d$participant == levels(d$participant)[1] & d$condition == "B"
  d$y[out] <- d$y[out] + 3
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d)
  loso <- diag_loso(m, "participant", terms = "conditionB")
  worst <- names(which.max(abs(loso$std_change[, "conditionB"])))
  expect_equal(worst, levels(d$participant)[1])
  expect_output(print(loso), "more than 1 SE")
  expect_error(diag_loso(m, "participant", terms = "nope"), "not")
})

test_that("diag_agreement detects a missing time course", {
  d <- sim_blocks("within", n_participants = 15, n_trials = 10, n_time = 20,
                  effect_condition = 0.5, seed = 3)
  m_ok <- lme4::lmer(y ~ condition + time + (1 | participant), data = d)
  ag_ok <- diag_agreement(m_ok, "condition", "participant")
  expect_false(ag_ok$flagged)
  expect_equal(ag_ok$estimate, ag_ok$mean_diff, tolerance = 0.05)
  expect_output(print(ag_ok), "agrees")

  # A covariate associated with condition makes the adjusted coefficient
  # differ from the raw participant-level differences
  set.seed(1)
  d$x <- ifelse(d$condition == "B", 1, 0) + stats::rnorm(nrow(d), 0, 0.1)
  d$y2 <- d$y + 2 * d$x
  m_adj <- lme4::lmer(y2 ~ condition + x + time + (1 | participant), data = d)
  ag_adj <- diag_agreement(m_adj, "condition", "participant")
  expect_true(ag_adj$flagged)
  expect_output(print(ag_adj), "differs")
  expect_equal(nrow(tidy(ag_adj)), 15)
})

test_that("diag_agreement validates its input", {
  d <- sim_blocks("within", "trial", n_participants = 4, n_trials = 4,
                  seed = 1)
  d$three <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  m <- lm(y ~ three, data = d)
  expect_error(diag_agreement(m, "three", "participant"), "two levels")
  expect_error(diag_agreement(m, "three", "site"), "not found")
})

test_that("influence and agreement plots are ggplot objects", {
  skip_if_not_installed("ggplot2")
  d <- sim_blocks("within", "trial", n_participants = 5, n_trials = 6,
                  seed = 4)
  m <- lm(y ~ condition + participant, data = d)
  expect_s3_class(plot(diag_loso(m, "participant", "conditionB")), "ggplot")
  expect_s3_class(plot(diag_agreement(m, "condition", "participant")),
                  "ggplot")
})
