# lmmresample

<!-- badges: start -->
[![R-CMD-check](https://github.com/ricardo-martinez-flores/lmmresample/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ricardo-martinez-flores/lmmresample/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Resampling inference for mixed models of densely sampled physiological time series.**

`lmmresample` provides permutation and bootstrap inference for linear and
generalized linear mixed-effects models (`lme4::lmer()`, `lme4::glmer()`)
fitted to densely sampled signals such as pupillometry or EEG, where
thousands of autocorrelated samples per trial make standard Wald or
Satterthwaite inference unreliable.

> **Status:** under active development. The interface is not yet stable.

## Planned for version 0.1.0

- **Permutation tests** that refit the full time-series model, under
  user-declared exchangeability: trials within participants, condition labels
  between participants, or custom strata.
- **Calibration of the permutation null** by simulating from the user's reduced
  model, optionally with AR(1) residuals, to verify type I error before
  interpreting a result.
- **Family-wise error control** across several model terms with the max-t
  procedure, using one shared permutation for all terms.
- **Bootstrap confidence intervals** with a choice of resampling scheme
  (cluster/case or wild) and interval type (percentile, BCa).
- **Diagnostics:** leave-one-subject-out influence and a coefficient agreement
  check that flags when the model coefficient diverges from the mean
  participant-level difference.
- Parallel computation through the `future` framework, with reproducible
  seeds.

Later versions: Freedman–Lane / ter Braak permutation with covariates, crossed
random-effects designs, F/likelihood-ratio statistics for multi-level factors,
and `glmmTMB` models.

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
