# lmmresample 0.0.0.9000

* `diag_acf()` now returns a bias-corrected AR(1) coefficient (`phi`).
  Centring each short series biased the mean autocorrelation strongly towards
  zero (e.g. 0.64 for a true 0.92 with 20 samples per series); the
  coefficient is now chosen so that the expected autocorrelations of centred
  AR(1) series of the observed lengths match the observed ones over all
  lags. The uncorrected lag-1 mean is kept as `phi_raw`, and
  `correct = FALSE` restores the previous behaviour. `perm_calibrate()`
  therefore simulates residuals with realistic dependence.

* Refits of mixed models start from the variance parameters of the original
  fit and skip lme4's derivative-based convergence check, roughly halving
  the time per refit on large data; `options(lmmresample.fast = FALSE)`
  restores the default settings.

* Fixes from an independent review:
  - Freedman-Lane: fitted values of the reduced model now keep only random
    effects grouped above the resampled units, so effects can no longer be
    absorbed by participant or trial intercepts when those units are
    permuted.
  - Freedman-Lane permutation of residuals between units now requires `time`
    and matches samples by time, refusing units with different time points.
  - `perm_maxt()` with min-p includes the observed data in the reference set
    of each test, which removes an anti-conservative bias for small `B`.
  - `data` supplied in a different order from the fitted model, and models
    with weights or offsets, are refused.
  - Series identifiers that repeat across participants are refused in
    `diag_acf()` and `perm_calibrate()`.
  - Default coefficients in `boot_ci()` and `diag_loso()` exclude terms of
    the cluster factor; `diag_agreement()` handles any contrast coding;
    groupings that include the cluster are no longer reported as crossed.
  - Random-effects terms left empty in the reduced model are dropped;
    calibration is reproducible with any RNG kind; the sharp null warns when
    there is no random slope to remove.

* New `boot_ci()` computes bootstrap confidence intervals for fixed effects by
  resampling whole clusters (case bootstrap, with repeated clusters and
  nested units relabelled) or by cluster wild bootstrap with Rademacher
  weights. Percentile, BCa (jackknife acceleration over clusters) and basic
  intervals; `confint()` returns any type without refitting. Random effects
  crossed with the cluster trigger a warning.

* `perm_maxt()` gains `method = "freedman-lane"`, so families of tests can
  include interactions and covariate-adjusted effects, and different terms
  in different models, and `combine = "min-p"`, which combines statistics on
  different scales (e.g. t and chi-square). The adjustment procedure is now
  chosen with `adjust`.

* New `diag_loso()` refits the model leaving out each participant in turn and
  flags influential participants; `diag_agreement()` compares the model
  coefficient of a two-level condition with the mean participant-level
  difference.

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
