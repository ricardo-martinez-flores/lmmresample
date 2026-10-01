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
`lm()`, `glm()`, `lme4::lmer()` and `lme4::glmer()`.

> **Status:** under active development. The interface is not yet stable.

## Planned for version 0.1.0

- **Permutation tests** that refit the full model, under
  user-declared exchangeability: trials within participants, condition labels
  between participants, or custom strata.
- **Calibration of the permutation null** by simulating from the user's reduced
  model, optionally with AR(1) residuals, to verify type I error before
  interpreting a result.
- **Family-wise error control** across several model terms with the max-t
  procedure, using one shared permutation for all terms, for example across
  several features extracted from the same trials.
- **Bootstrap confidence intervals** with a choice of resampling scheme
  (cluster/case or wild) and interval type (percentile, BCa).
- **Diagnostics:** leave-one-subject-out influence and a coefficient agreement
  check that flags when the model coefficient diverges from the mean
  participant-level difference.
- **Diagnostic plots** for every result: null and bootstrap distributions,
  calibration of p-values, leave-one-subject-out influence and residual
  autocorrelation.
- Parallel computation through the `future` framework, with reproducible
  seeds.

Planned for version 0.2: Freedman–Lane / ter Braak permutation with
covariates (including tests of interactions), F / likelihood-ratio statistics
for multi-level factors, and `glmmTMB` models.

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
