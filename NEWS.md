# lmmresample 0.0.0.9000

* New diagnostics: `diag_acf()` estimates the autocorrelation of residuals
  within series (its lag-1 value can be passed to `perm_calibrate()`), and
  `diag_timecourse()` shows mean residuals over time by condition to detect a
  misspecified time course or a missing condition by time interaction.

* The package now supports linear models only (`lm()` and `lme4::lmer()`).
  Generalized linear (mixed) models are refused with an informative message,
  because residual permutation does not apply to them.

* `perm_test()` and `perm_calibrate()` gain `method = "freedman-lane"`:
  residuals of the reduced model are permuted between whole units or
  sign-flipped by block, which allows tests of interactions (e.g.
  `condition:time`), of terms with several coefficients (joint Wald
  chi-square) and of effects in the presence of nuisance terms. Random slopes
  of the tested term are removed from the reduced model. `method = "auto"`
  chooses relabelling for simple main effects and Freedman-Lane otherwise.

* `sim_blocks()` gains `effect_onset` (condition effects that appear later in
  the trial) and `ar1_trials` (correlation between consecutive trials).

* New `perm_maxt()` and `perm_spec()` test the same variable in several models
  at once (e.g. one model per feature) and control the family-wise error rate
  with the Westfall-Young max-t procedure (single-step or step-down), using
  one shared relabelling of units per permutation.

* New `perm_calibrate()` checks the type I error of a permutation test for the
  user's own design by simulating data under the null from the fitted model
  (mean or sharp null, optional AR(1) residuals within series), and compares
  it with the model's Wald test. Includes p-value calibration and
  rejection-rate plots.

* New `perm_test()` tests a fixed-effect term by permuting it across whole
  units (trials or participants) and refitting the full model. Supports
  `lm()` and `lme4::lmer()` models, runs in parallel
  through the future framework, and records failed, non-convergent and
  singular refits. Results have `print()`, `summary()`, `tidy()` and `plot()`
  methods, including a trace plot of the running p-value.

* New `exch_within()`, `exch_signflip()` and `exch_free()` declare which units
  are exchangeable: within blocks (e.g. trials within participants), by
  swapping the two levels of the tested variable for whole blocks
  (sign-flipping, for the mean effect when participants differ in their
  effect), freely, or within strata.

* `sim_blocks()` gains `prop_condition` (unbalanced conditions, e.g. oddball
  designs), `residual_df` (heavy-tailed residuals), and per-participant trial
  counts through a vector `n_trials`.

* New `sim_blocks()` simulates block-structured repeated measurements (full
  time series, trial-level or participant-level data) for within-participant,
  between-group and mixed designs, with known effects, random intercepts and
  slopes, and AR(1) residuals within trials.
