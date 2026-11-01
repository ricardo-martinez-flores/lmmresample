# lmmresample

<!-- badges: start -->
[![R-CMD-check](https://github.com/ricardo-martinez-flores/lmmresample/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ricardo-martinez-flores/lmmresample/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Resampling inference for mixed models of densely sampled physiological time series.**

`lmmresample` provides permutation and bootstrap inference for models fitted
to repeated measurements organised in blocks (trials, epochs or
participants). Resampling always operates on whole blocks and refits a single
model on the full data, so dependence within each block, such as the strong
autocorrelation of samples within a trial, is carried into the reference
distribution without having to be modelled. Where hundreds of autocorrelated
samples per trial make standard Wald or Satterthwaite inference unreliable,
this gives valid tests and intervals.

It applies to any block-structured signal, including pupillometry, eye
movements, EEG, fNIRS, heart rate and HRV, electrodermal activity, EMG and
kinematics, and at any level of aggregation: the full time series, features
extracted per trial, or one value per participant. Supported models are
`lm()` and `lme4::lmer()`.

> **Status:** under active development. The interface is not yet stable.

## Features

- **Permutation tests** that refit the full model, under user-declared
  exchangeability: trials within participants, condition labels between
  participants, sign-flipping of whole participants, or custom strata.
- **Freedman–Lane permutation** of the residuals of the reduced model, for
  tests of interactions (e.g. condition × time), of terms with several
  coefficients and of effects adjusted for covariates.
- **Calibration of the permutation test** by simulating from the user's own
  model under the null hypothesis, optionally with AR(1) residuals, to check
  type I error for the actual design before interpreting a result.
- **Family-wise error control** across several outcomes with the max-t
  procedure, using one shared permutation for all tests, for example across
  several features extracted from the same trials.
- **Pairwise comparisons** between the levels of a factor (post hoc tests),
  in level or in time course, with Freedman–Lane tests that keep the other
  levels in the model, family-wise adjustment across pairs and optional
  bootstrap intervals of each difference.
- **Bootstrap confidence intervals** by resampling whole participants (case
  or wild bootstrap), with expanded percentile (the default, corrected for
  few participants), percentile, BCa and basic intervals.
- **Diagnostics:** residual autocorrelation, residual time courses by
  condition, leave-one-participant-out influence, and agreement between the
  model coefficient and participant-level differences.
- **Simulation of block-structured data** with known effects, for examples,
  power analysis and checking an analysis pipeline.
- **Plots** for every result and parallel computation through the `future`
  framework, with reproducible seeds.

## Installation

Development version from GitHub:

```r
# install.packages("pak")
pak::pak("ricardo-martinez-flores/lmmresample")
```

## Example

A permutation test of a condition × time interaction in a mixed model of
full time series, and bootstrap intervals for its coefficients:

```r
library(lmmresample)
library(lme4)
library(splines)

d <- sim_blocks("within", level = "timeseries", n_participants = 20,
                n_trials = 20, n_time = 50, sampling_rate = 25,
                effect_condition = 0.6, effect_onset = 0.5, seed = 2)
m <- lmer(y ~ condition * ns(time, 4) + (1 | participant), data = d)

# Whole trials are permuted; samples within a trial are matched by time
perm_test(m, "condition:ns(time, 4)", unit = "trial",
          exchange = exch_signflip("participant"), time = "time")

# Whole participants are resampled
boot_ci(m, cluster = "participant",
        terms = grep("^conditionB:", names(fixef(m)), value = TRUE))

# Does the test hold its error rate for this design?
ac <- diag_acf(m, series = "trial", time = "time")
perm_calibrate(m, "condition:ns(time, 4)", unit = "trial",
               exchange = exch_signflip("participant"), ar1 = ac,
               series = "trial", time = "time")
```

Computation runs in parallel with `future::plan("multisession")`.

## Learning more

- `vignette("lmmresample")`: getting started with a trial-level analysis.
- `vignette("time-courses")`: full time series, diagnostics and calibration.
- `vignette("multiple-comparisons")`: families of tests and pairwise
  comparisons.

## Authors

- Ricardo Martínez-Flores (maintainer) — [ORCID 0000-0002-8435-2710](https://orcid.org/0000-0002-8435-2710)
- Hans Supèr — [ORCID 0000-0001-9328-0096](https://orcid.org/0000-0001-9328-0096)

## License

MIT © lmmresample authors
