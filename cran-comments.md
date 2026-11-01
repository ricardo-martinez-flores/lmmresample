## Test environments

* macOS 26.5 (aarch64), R 4.4.1
* Ubuntu 22.04, R 4.6.1
* win-builder, R-devel (2026-10-31 r90820 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION: electrodermal,
  electromyography, exchangeability, fNIRS, pupillometry. These are
  technical terms spelled correctly (fNIRS: functional near-infrared
  spectroscopy).

## Vignettes

The vignettes are precomputed: they run permutation and bootstrap analyses
that take several minutes, so the shipped .Rmd files contain the results
obtained with the code shown, and building them does not repeat the
computation.
