test_that("constructors validate their arguments", {
  expect_s3_class(exch_within("participant"), "lmmr_exchange")
  expect_s3_class(exch_free(), "lmmr_exchange")
  expect_s3_class(exch_free(strata = "site"), "lmmr_exchange")
  expect_error(exch_within(), "block")
  expect_error(exch_within(c("a", "b")), "single column name")
  expect_error(exch_free(strata = 1), "single column name")
})

test_that("exchangeability is printed in words", {
  expect_output(print(exch_within("participant")), "within `participant`")
  expect_output(print(exch_free()), "freely")
  expect_output(print(exch_free("site")), "strata `site`")
})

test_that("permutations stay within blocks and preserve labels", {
  block <- rep(c("a", "b", "c"), times = c(4, 6, 1))
  values <- c("x", "y", "x", "y", "x", "x", "y", "y", "x", "y", "x")
  set.seed(1)
  perms <- generate_permutations(block, B = 50)

  expect_equal(dim(perms), c(length(block), 50))
  for (b in seq_len(50)) {
    expect_identical(block[perms[, b]], block)
    for (blk in unique(block)) {
      in_blk <- block == blk
      expect_identical(sort(values[perms[in_blk, b]]), sort(values[in_blk]))
    }
  }
  expect_true(all(perms[block == "c", ] == which(block == "c")))
})

test_that("unit table validates the design", {
  d <- sim_blocks("within", "trial", n_participants = 4, n_trials = 4,
                  seed = 1)
  ub <- build_units(d, "condition", "trial", exch_within("participant"))
  expect_equal(nrow(ub$units), nrow(d))
  expect_equal(length(unique(ub$units$block)), 4)

  expect_error(
    build_units(d, "condition", "participant", exch_free()),
    "constant within each"
  )
  expect_error(
    build_units(d, "condition", "trial", exch_within("site")),
    "not found"
  )
  expect_error(
    build_units(d, "condition", "trial", "within"),
    "exch_within"
  )

  d$site <- rep(c("s1", "s2"), length.out = nrow(d))
  expect_error(
    build_units(d, "condition", "participant",
                exch_free(strata = "site")),
    "constant within each"
  )
  db <- sim_blocks("between", "trial", n_participants = 4, n_trials = 4,
                   seed = 1)
  db$site <- rep(c("s1", "s2"), length.out = nrow(db))
  expect_error(
    build_units(db, "group", "participant", exch_free(strata = "site")),
    "single stratum"
  )
})

test_that("the number of distinct permutations is counted", {
  n <- count_permutations(c("x", "x", "y", "y"), rep("all", 4))
  expect_equal(exp(n$log_n), choose(4, 2))
  n <- count_permutations(c("x", "y", "x", "y"), c("a", "a", "b", "b"))
  expect_equal(exp(n$log_n), 4)
  n <- count_permutations(c("x", "x", "x", "y"), c("a", "a", "b", "b"))
  expect_equal(n$uninformative, 1)
})

test_that("sign-flipping swaps the two values for whole blocks", {
  codes <- c(1L, 2L, 1L, 2L, 2L, 1L)
  block <- c("a", "a", "b", "b", "c", "c")
  set.seed(2)
  relabels <- generate_relabels(codes, block, "signflip", 200)
  for (b in seq_len(200)) {
    for (blk in unique(block)) {
      in_blk <- block == blk
      same <- all(relabels[in_blk, b] == codes[in_blk])
      flipped <- all(relabels[in_blk, b] == 3L - codes[in_blk])
      expect_true(same || flipped)
    }
  }
  flipped_a <- mean(relabels[1, ] != codes[1])
  expect_gt(flipped_a, 0.35)
  expect_lt(flipped_a, 0.65)
})

test_that("sign-flipping is described and validated", {
  expect_output(print(exch_signflip("participant")), "sign-flipped")
  expect_error(exch_signflip(), "block")
  ex <- exch_signflip("participant")
  expect_error(encode_values(c("a", "b", "c"), ex), "exactly two values")
  enc <- encode_values(factor(c("B", "A", "B"), levels = c("A", "B")), ex)
  expect_equal(enc$codes, c(2L, 1L, 2L))
})

test_that("sign-flipping counts 2^(blocks - 1) distinct relabellings", {
  n <- count_permutations(rep(c("x", "y"), 5), rep(letters[1:5], each = 2),
                          "signflip")
  expect_equal(exp(n$log_n), 2^4)
})
