# lmmresample 0.0.0.9000

* New `perm_calibrate()` checks the type I error of a permutation test for the
  user's own design by simulating data under the null from the fitted model
  (mean or sharp null, optional AR(1) residuals within series), and compares
  it with the model's Wald test. Includes p-value calibration and
  rejection-rate plots.

* New `perm_test()` tests a fixed-effect term by permuting it across whole
  units (trials or participants) and refitting the full model. Supports
  `lm()`, `glm()`, `lme4::lmer()` and `lme4::glmer()` models, runs in parallel
  through the future framework, and records failed, non-convergent and
  singular refits. Results have `print()`, `summary()`, `tidy()` and `plot()`
  methods, including a trace plot of the running p-value.

* New `exch_within()`, `exch_signflip()` and `exch_free()` declare which units
  are exchangeable: within blocks (e.g. trials within participants), by
  swapping the two levels of the tested variable for whole blocks
  (sign-flipping, for the mean effect when participants differ in their
  effect), freely, or within strata.

* `sim_blocks()` gains `prop_condition` (unbalanced conditions, e.g. oddball
  designs) and `residual_df` (heavy-tailed residuals).

* New `sim_blocks()` simulates block-structured repeated measurements (full
  time series, trial-level or participant-level data) for within-participant,
  between-group and mixed designs, with known effects, random intercepts and
  slopes, and AR(1) residuals within trials.
