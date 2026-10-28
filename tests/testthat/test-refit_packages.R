test_that("refit_packages includes the fitting package and formula functions", {
  skip_if_not_installed("lme4")
  d <- data.frame(y = rnorm(40), x = seq_len(40) / 40,
                  g = factor(rep(1:8, each = 5)))
  ns <- splines::ns
  m <- lme4::lmer(y ~ ns(x, 2) + (1 | g), data = d)
  pk <- refit_packages(m)
  expect_true("lme4" %in% pk)
  expect_true("splines" %in% pk)
  expect_false("lmmresample" %in% pk)
})

test_that("rejection_table handles simulations without valid p-values", {
  tab <- rejection_table(list(Permutation = c(NA_real_, NA_real_),
                              Wald = c(0.01, 0.5)), 0.05)
  expect_true(is.na(tab$rejection[1]))
  expect_equal(tab$n, c(0L, 2L))
  expect_equal(tab$rejection[2], 0.5)
})
