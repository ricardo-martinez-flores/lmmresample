boot_data <- function(seed = 1, ...) {
  sim_blocks("within", "trial", n_participants = 10, n_trials = 8,
             effect_condition = 0.5, seed = seed, ...)
}

test_that("boot_ci returns a complete lmmr_boot object", {
  d <- boot_data()
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d)
  bt <- boot_ci(m, "participant", ci = "bca", B = 49, seed = 1)
  expect_s3_class(bt, "lmmr_boot")
  expect_equal(colnames(bt$boot), "conditionB")
  expect_equal(nrow(bt$boot), 49)
  expect_equal(nrow(bt$jackknife), 10)
  iv <- bt$intervals
  expect_true(iv$conf.low < iv$estimate && iv$estimate < iv$conf.high)
  expect_output(print(bt), "BCa intervals")
  expect_output(print(bt), "exploratory")
})

test_that("intervals follow their definitions", {
  d <- boot_data()
  m <- lm(y ~ condition + participant, data = d)
  bt <- boot_ci(m, "participant", terms = "conditionB", ci = "bca", B = 99,
                seed = 2)
  b <- bt$boot[, 1]
  est <- bt$estimate[[1]]
  q <- stats::quantile(b, c(0.025, 0.975), names = FALSE, type = 6)
  expect_equal(unname(confint(bt, type = "percentile")[1, ]), q)
  expect_equal(unname(confint(bt, type = "basic")[1, ]), 2 * est - rev(q))

  jack <- bt$jackknife[, 1]
  z0 <- stats::qnorm(mean(b < est) + 0.5 * mean(b == est))
  dd <- mean(jack) - jack
  acc <- sum(dd^3) / (6 * sum(dd^2)^1.5)
  z <- stats::qnorm(c(0.025, 0.975))
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - acc * (z0 + z)))
  expect_equal(unname(confint(bt, type = "bca")[1, ]),
               stats::quantile(b, adj, names = FALSE, type = 6))
  expect_equal(colnames(confint(bt, level = 0.9)), c("5 %", "95 %"))
})

test_that("case resampling relabels repeated clusters and nested units", {
  d <- sim_blocks("within", n_participants = 6, n_trials = 4, n_time = 4,
                  seed = 3)
  m <- suppressMessages(
    lme4::lmer(y ~ condition + (1 | participant) + (1 | trial), data = d)
  )
  nested <- nested_grouping(m, d, "participant")
  expect_equal(nested, "trial")
  rs <- case_resampler(d, "participant", nested)
  draws <- matrix(c(1, 1, 2, 3, 4, 5), ncol = 1)
  new <- rs$build(draws, 1)
  expect_equal(nlevels(new$participant), 6)
  expect_equal(nlevels(new$trial), 6 * 4)
  expect_equal(nrow(new), nrow(d))
})

test_that("crossed grouping factors trigger a warning", {
  d <- sim_blocks("within", "trial", n_participants = 6, n_trials = 6,
                  seed = 4)
  d$item <- factor(rep(paste0("i", 1:6), length.out = nrow(d)))
  m <- suppressMessages(
    lme4::lmer(y ~ condition + (1 | participant) + (1 | item), data = d)
  )
  expect_warning(boot_ci(m, "participant", B = 9, seed = 1), "crossed")
})

test_that("the wild bootstrap keeps the design and resamples signs", {
  d <- boot_data()
  m <- lme4::lmer(y ~ condition + (1 | participant), data = d)
  rs <- wild_resampler(m, d, "participant")
  set.seed(1)
  draws <- rs$draw(3)
  new <- rs$build(draws, 1)
  expect_identical(new$condition, d$condition)
  xb <- stats::predict(m, re.form = NA)
  e <- d$y - xb
  ratio <- (new$y - xb) / e
  per_cluster <- tapply(round(ratio, 8), d$participant,
                        function(r) length(unique(r)))
  expect_true(all(per_cluster == 1))
  expect_true(all(abs(round(ratio, 8)) == 1))

  bt <- boot_ci(m, "participant", resample = "wild", B = 29, seed = 1)
  expect_equal(bt$resample, "wild")
})

test_that("results are reproducible and inputs validated", {
  d <- boot_data()
  m <- lm(y ~ condition + participant, data = d)
  b1 <- boot_ci(m, "participant", terms = "conditionB", B = 19, seed = 5)
  b2 <- boot_ci(m, "participant", terms = "conditionB", B = 19, seed = 5)
  expect_identical(b1$boot, b2$boot)
  expect_error(boot_ci(m, "site"), "not found")
  expect_error(boot_ci(m, "participant", terms = "nope"), "not")
  expect_error(boot_ci(m, "participant", level = 1), "smaller than 1")
  d2 <- d[d$participant %in% levels(d$participant)[1:2], ]
  m2 <- lm(y ~ condition, data = droplevels(d2))
  expect_error(boot_ci(m2, "participant", B = 9), "At least 3")
})

test_that("tidy and plot methods work", {
  d <- boot_data()
  m <- lm(y ~ condition + participant, data = d)
  bt <- boot_ci(m, "participant", terms = "conditionB", B = 19, seed = 1)
  expect_equal(nrow(tidy(bt)), 1)
  expect_equal(nrow(tidy(bt, type = "replicates")), 19)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(bt), "ggplot")
})

test_that("bootstrap intervals cover the true effect", {
  skip_on_cran()
  n_sim <- 40
  cover <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    d <- sim_blocks("within", "trial", n_participants = 20, n_trials = 10,
                    effect_condition = 0.5, sd_slope = 0.3, seed = i)
    m <- suppressMessages(
      lme4::lmer(y ~ condition + (1 + condition | participant), data = d)
    )
    ci <- confint(suppressWarnings(
      boot_ci(m, "participant", terms = "conditionB", B = 199, seed = i)
    ))
    cover[i] <- ci[1, 1] <= 0.5 && 0.5 <= ci[1, 2]
  }
  expect_gt(mean(cover), 0.8)
})

test_that("BCa without the jackknife gives an informative error", {
  d <- boot_data()
  m <- lm(y ~ condition + participant, data = d)
  bt <- boot_ci(m, "participant", B = 19, seed = 1)
  expect_null(bt$jackknife)
  expect_error(confint(bt, type = "bca"), "ci = \"bca\"")
})

test_that("expanded percentile is the default and widens with few clusters", {
  d <- sim_blocks("within", "trial", n_participants = 8, n_trials = 10,
                  effect_condition = 0.5, sd_slope = 0.3, seed = 2)
  m <- suppressMessages(
    lme4::lmer(y ~ condition + (1 + condition | participant), data = d)
  )
  bt <- suppressWarnings(boot_ci(m, "participant", terms = "conditionB",
                                 B = 199, seed = 1))
  expect_identical(bt$ci, "expanded")
  pe <- confint(bt, type = "percentile")
  ex <- confint(bt)
  expect_lt(ex[1, 1], pe[1, 1])
  expect_gt(ex[1, 2], pe[1, 2])
  b <- bt$boot[, "conditionB"]
  G <- 8
  z <- sqrt(G / (G - 1)) * stats::qt(0.975, G - 1)
  q <- stats::quantile(b[!is.na(b)], stats::pnorm(c(-z, z)), names = FALSE,
                       type = 6)
  expect_equal(unname(ex[1, ]), q)
})

strata_data <- function(n_control = 8, n_patient = 4, seed = 3) {
  d <- sim_blocks("within", "trial", n_participants = n_control + n_patient,
                  n_trials = 12, effect_condition = 0.5, sd_slope = 0.2,
                  seed = seed)
  parts <- levels(d$participant)
  grp <- rep(c("control", "patient"), c(n_control, n_patient))
  d$group <- factor(grp[match(d$participant, parts)])
  d
}

test_that("stratified case resampling keeps the clusters per stratum", {
  d <- strata_data()
  m <- suppressMessages(
    lme4::lmer(y ~ group * condition + (1 | participant), data = d)
  )
  st <- cluster_strata(d, "participant", "group")
  rs <- case_resampler(d, "participant", character(0), st)
  draws <- with_seed(1, rs$draw(200))
  ids <- unique(as.character(d$participant))
  # every position is filled with a cluster of its own stratum
  expect_true(all(st[draws] == rep(st, times = 200)))
  nd <- rs$build(draws, 1)
  per_group <- tapply(as.character(nd$participant), nd$group,
                      function(x) length(unique(x)))
  expect_equal(as.numeric(per_group[c("control", "patient")]), c(8, 4))
  bt <- suppressWarnings(boot_ci(m, "participant", strata = "group",
                                 terms = "grouppatient:conditionB",
                                 B = 49, seed = 1))
  expect_equal(bt$n_strata, 2L)
  expect_identical(bt$strata, "group")
  expect_equal(bt$B_used, 49)
  b <- bt$boot[, 1]
  G <- 12; S <- 2
  z <- sqrt(G / (G - S)) * stats::qt(0.975, G - S)
  q <- stats::quantile(b, stats::pnorm(c(-z, z)), names = FALSE, type = 6)
  expect_equal(unname(confint(bt)[1, ]), q)
  expect_output(print(bt), "within 2 strata")
})

test_that("strata must be constant within clusters and not too small", {
  d <- strata_data()
  m <- suppressMessages(
    lme4::lmer(y ~ group * condition + (1 | participant), data = d)
  )
  expect_error(boot_ci(m, "participant", strata = "condition", B = 9),
               "constant within each")
  expect_error(boot_ci(m, "participant", strata = "nothere", B = 9),
               "not found")
  d2 <- strata_data(n_control = 10, n_patient = 1)
  m2 <- suppressMessages(
    lme4::lmer(y ~ group * condition + (1 | participant), data = d2)
  )
  expect_error(boot_ci(m2, "participant", strata = "group", B = 9),
               "at least 2 clusters")
})

test_that("strata are ignored by the wild bootstrap", {
  d <- strata_data()
  m <- suppressMessages(
    lme4::lmer(y ~ group * condition + (1 | participant), data = d)
  )
  a <- suppressWarnings(boot_ci(m, "participant", resample = "wild",
                                strata = "group", B = 19, seed = 2))
  b <- suppressWarnings(boot_ci(m, "participant", resample = "wild",
                                B = 19, seed = 2))
  expect_equal(a$boot, b$boot)
  expect_null(a$strata)
  expect_equal(a$n_strata, 1L)
})
