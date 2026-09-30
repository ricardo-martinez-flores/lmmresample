#' lmmresample: Resampling Inference for Mixed Models of Densely Sampled Time Series
#'
#' Permutation and bootstrap inference for mixed-effects models fitted to
#' densely sampled physiological time series (pupillometry, EEG).
#'
#' The package is under active development. Planned for version 0.1.0:
#'
#' * Permutation tests under user-declared exchangeability (trials within
#'   participants, labels between participants, or custom strata).
#' * Simulation-based calibration of the permutation null.
#' * Family-wise error control across several model terms (max-t).
#' * Cluster bootstrap confidence intervals (percentile and BCa), with wild
#'   bootstrap as an alternative resampling scheme.
#' * Leave-one-subject-out diagnostics and a coefficient agreement check.
#'
#' @keywords internal
"_PACKAGE"
