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

Available now:

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
- **Simulation of block-structured data** with known effects, for examples,
  power analysis and checking an analysis pipeline.
- **Plots** for every result and parallel computation through the `future`
  framework, with reproducible seeds.

In development:

- **Bootstrap confidence intervals** with a choice of resampling scheme
  (cluster or wild) and interval type (percentile, BCa).
- **Diagnostics:** leave-one-subject-out influence, residual autocorrelation,
  residual time courses by condition, and a coefficient agreement check.

## Installation

Development version from GitHub:

```r
# install.packages("pak")
pak::pak("ricardo-martinez-flores/lmmresample")
```

## Authors

- Ricardo Martínez-Flores (maintainer) — [ORCID 0000-0002-8435-2710](https://orcid.org/0000-0002-8435-2710)
- Hans Supèr

## License

MIT © lmmresample authors
